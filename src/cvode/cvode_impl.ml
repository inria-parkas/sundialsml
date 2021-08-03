(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

open Sundials

(* Hack to ensure that Sundials.c_init_module is executed so that the global
   exceptions are properly registered. *)
let e = Sundials.RecoverableFailure

(* Types shared between Cvode, Cvodes, Cvode_bbd, and Cvodes_bbd.  *)

(* This module's purpose is just to define types shared between Cvode,
   Cvodes, Cvode_bbd, and Cvode_bbds.  However, the session and
   session_linear_solver types must have their implementation exposed
   to all four modules yet be abstract to code that lives outside of
   the sundialsml library.

   To satisfy this requirement, we define the type in a separate
   module - this module - that exports everything, then omit
   Cvode_impl.cmi during installation.  This way, the compiled library
   doesn't divulge implementation details.

   Unfortunately, session depends on 90% of all data types used in
   those modules, so we are forced to declare almost all types here.
   To avoid repeating the types' definitions in Cvode, Cvodes,
   Cvode_bbd, and Cvodes_bbd, we "include Cvode_impl" from each
   module.  However, some submodules' types clash (for example
   Cvode.Spils.prec_solve_fn and Cvodes.Adjoint.Spils.prec_solve_fn),
   so we need to reproduce to some extent the submodules present in
   Cvode and Cvodes.  Hence the behemoth you see below.

 *)

module LSI = Sundials_LinearSolver_impl
module NLSI = Sundials_NonlinearSolver_impl

(*
 * NB: The order of variant constructors and record fields is important!
 *     If these types are changed or augmented, the corresponding declarations
 *     in cvode_ml.h (and code in cvode_ml.c) must also be updated.
 *)

type ('data, 'kind) nvector = ('data, 'kind) Nvector.t

type 'a double = 'a * 'a
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

  (* These fields are accessed from cvode_ml.c *)
  type 'm jac_callback =
    {
      jacfn: 'm jac_fn;
      mutable jmat : 'm option (* Not used in Sundials >= 3.0.0 *)
    }

  type 'm linsys_fn =
    (RealArray.t triple, RealArray.t) jacobian_arg -> 'm -> bool -> float -> bool

  let no_callback = fun _ _ -> Sundials_impl.crash "no direct callback"

  let no_linsysfn = fun _ _ _ _ ->
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

module CvodeBbdParamTypes = struct
  type 'a local_fn = float -> 'a -> 'a -> unit
  type 'a comm_fn = float -> 'a -> unit
  type 'a precfns =
    {
      local_fn : 'a local_fn;
      comm_fn  : 'a comm_fn option;
    }
end

module CvodeBbdTypes = struct
  type bandwidths =
    {
      mudq    : int;
      mldq    : int;
      mukeep  : int;
      mlkeep  : int;
    }
end

(* Sensitivity *)

module QuadratureTypes = struct
  type 'a quadrhsfn = float -> 'a -> 'a -> unit
end

module SensitivityTypes = struct
  type 'd sensrhsfn_args =
    {
      t : float;
      y : 'd;
      y' : 'd;
      tmp : 'd double;
    }

  type 'a sensrhsfn_all =
    'a sensrhsfn_args
    -> 'a array
    -> 'a array
    -> unit

  type 'a sensrhsfn1 =
    int
    -> 'a sensrhsfn_args
    -> 'a
    -> 'a
    -> unit

  type 'a sensrhsfn =
      AllAtOnce of 'a sensrhsfn_all option
    | OneByOne of 'a sensrhsfn1 option

  module QuadratureTypes = struct
    type 'd quadsensrhsfn_args =
      {
        t : float;
        y : 'd;
        s : 'd array;
        yq' : 'd;
        tmp : 'd double;
      }

    type 'a quadsensrhsfn = 'a quadsensrhsfn_args -> 'a array -> unit
  end
end

module AdjointTypes' = struct
  type 'd brhsfn_args =
    {
      t : float;
      y : 'd;
      yb : 'd;
    }
  type 'a brhsfn_no_sens = 'a brhsfn_args -> 'a -> unit
  type 'a brhsfn_with_sens = 'a brhsfn_args -> 'a array -> 'a -> unit

  type 'a brhsfn =
      NoSens of 'a brhsfn_no_sens
    | WithSens of 'a brhsfn_with_sens

  module QuadratureTypes = struct
    type 'd bquadrhsfn_args =
      {
        t : float;
        y : 'd;
        yb : 'd;
      }
    type 'a bquadrhsfn_no_sens = 'a bquadrhsfn_args -> 'a -> unit
    type 'a bquadrhsfn_with_sens = 'a bquadrhsfn_args -> 'a array -> 'a -> unit
    type 'a bquadrhsfn =
        NoSens of 'a bquadrhsfn_no_sens
      | WithSens of 'a bquadrhsfn_with_sens
  end

  type ('t, 'a) jacobian_arg =
    {
      jac_t   : float;
      jac_y   : 'a;
      jac_yb  : 'a;
      jac_fyb : 'a;
      jac_tmp : 't
    }

  (* This is NOT the same as DirectTypes defined above.  This version
     refers to a different jacobian_arg, the one that was just
     defined.  *)
  module DirectTypes = struct
    type 'm jac_fn_no_sens
      = (RealArray.t triple, RealArray.t) jacobian_arg
        -> 'm
        -> unit

    type 'm jac_fn_with_sens
      = (RealArray.t triple, RealArray.t) jacobian_arg
        -> RealArray.t array
        -> 'm
        -> unit

    type 'm jac_fn =
      NoSens of 'm jac_fn_no_sens
    | WithSens of 'm jac_fn_with_sens

    (* These fields are accessed from cvodes_ml.c *)
    type 'm jac_callback_no_sens =
      {
        jacfn: 'm jac_fn_no_sens;
        mutable jmat : 'm option
      }

    let no_callback = fun _ _ -> Sundials_impl.crash "no direct callback"

    (* These fields are accessed from cvodes_ml.c *)
    type 'm jac_callback_with_sens =
      {
        jacfn_sens : 'm jac_fn_with_sens;
        mutable jmat : 'm option
      }

    type 'm linsys_fn_no_sens =
      (RealArray.t triple, RealArray.t) jacobian_arg
      -> 'm
      -> bool
      -> float
      -> bool

    type 'm linsys_fn_with_sens =
      (RealArray.t triple, RealArray.t) jacobian_arg
      -> RealArray.t array
      -> 'm
      -> bool
      -> float
      -> bool

    type 'm linsys_fn =
        LNoSens of 'm linsys_fn_no_sens
      | LWithSens of 'm linsys_fn_with_sens

    let no_linsysfn =
      LNoSens (fun _ _ _ _ ->
        Sundials_impl.crash "no linear system function no sens callback")

  end

  (* Ditto. *)
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

    type 'a precfns_no_sens =
      {
        prec_solve_fn : 'a prec_solve_fn;
        prec_setup_fn : 'a prec_setup_fn option;
      }

    type 'd jac_times_setup_fn_no_sens = (unit, 'd) jacobian_arg -> unit

    type 'a jac_times_vec_fn_no_sens =
      ('a, 'a) jacobian_arg
      -> 'a (* v *)
      -> 'a (* Jv *)
      -> unit

    (* versions with forward sensitivities *)

    type 'd prec_solve_fn_with_sens =
      (unit, 'd) jacobian_arg
      -> 'd prec_solve_arg
      -> 'd array
      -> 'd
      -> unit

    type 'd prec_setup_fn_with_sens =
      (unit, 'd) jacobian_arg
      -> 'd array
      -> bool
      -> float
      -> bool

    type 'a precfns_with_sens =
      {
        prec_solve_fn_sens : 'a prec_solve_fn_with_sens;
        prec_setup_fn_sens : 'a prec_setup_fn_with_sens option;
      }

    type 'd jac_times_setup_fn_with_sens =
      (unit, 'd) jacobian_arg -> 'd array -> unit

    type 'd jac_times_vec_fn_with_sens =
      ('d, 'd) jacobian_arg
      -> 'd array
      -> 'd
      -> 'd
      -> unit
  end
end

module CvodesBbdParamTypes = struct
  type 'a local_fn = 'a AdjointTypes'.brhsfn_args -> 'a -> unit
  type 'a comm_fn = 'a AdjointTypes'.brhsfn_args -> unit
  type 'a precfns =
    {
      local_fn : 'a local_fn;
      comm_fn  : 'a comm_fn option;
    }
end
module CvodesBbdTypes = CvodeBbdTypes

type cvode_mem
type c_weak_ref

type 'a rhsfn = float -> 'a -> 'a -> unit
type 'a rootsfn = float -> 'a -> RealArray.t -> unit
type error_handler = Util.error_details -> unit
type 'a error_weight_fun = 'a -> 'a -> unit
type 'd proj_fn = float -> 'd -> 'd -> float -> 'd -> unit

(* Session: here comes the big blob.  These mutually recursive types
   cannot be handed out separately to modules without menial
   repetition, so we'll just have them all here, at the top of the
   Types module.

   The ls_solver field only exists to ensure that the linear solver is not
   garbage collected while still being used by a session.
*)

(* Fields must be given in the same order as in cvode_session_index *)
type ('a, 'kind) session = {
  cvode      : cvode_mem;
  backref    : c_weak_ref;
  nroots     : int;
  checkvec   : (('a, 'kind) Nvector.t -> unit);

  mutable exn_temp     : exn option;

  rhsfn                : 'a rhsfn;
  mutable rootsfn      : 'a rootsfn;
  mutable errh         : error_handler;
  mutable errw         : 'a error_weight_fun;

  mutable projfn       : 'a proj_fn;
  mutable monitorfn    : ('a, 'kind) session -> unit;

  mutable ls_solver    : LSI.held_linear_solver;
  mutable ls_callbacks : ('a, 'kind) linsolv_callbacks;
  mutable ls_precfns   : 'a linsolv_precfns;

  mutable nls_solver   : ('a, 'kind, (('a, 'kind) session) NLSI.integrator)
                           NLSI.nonlinear_solver option;

  mutable sensext      : ('a, 'kind) sensext (* Used by Cvodes *)
}

(* Note: When compatibility with Sundials < 3.0.0 is no longer required,
         this type can be greatly simplified since we would no longer
         need to distinguish between different "direct" linear solvers.

   Note: The first field must always hold the Jacobian callback closure
         (it is accessed as Field(cb, 0) from cvode_ml.c.
         The second field must always hold the linear system callback closure
         (it is accessed as Field(cb, 1) from cvode_ml.c.
         This organization is used, despite the fact that the two callbacks
         are never used simultaneously, for simplicity.
*)
and ('a, 'kind) linsolv_callbacks =
  | NoCallbacks

  (* Diagonal *)
  | DiagNoCallbacks

  (* Dls *)
  | DlsDenseCallback of Matrix.Dense.t DirectTypes.jac_callback
                        * Matrix.Dense.t DirectTypes.linsys_fn
  | DlsBandCallback of Matrix.Band.t  DirectTypes.jac_callback
                        * Matrix.Band.t DirectTypes.linsys_fn

  | BDlsDenseCallback
      of Matrix.Dense.t AdjointTypes'.DirectTypes.jac_callback_no_sens
         * Matrix.Dense.t AdjointTypes'.DirectTypes.linsys_fn
  | BDlsDenseCallbackSens
      of Matrix.Dense.t AdjointTypes'.DirectTypes.jac_callback_with_sens
         * Matrix.Dense.t AdjointTypes'.DirectTypes.linsys_fn
  | BDlsBandCallback
      of Matrix.Band.t AdjointTypes'.DirectTypes.jac_callback_no_sens
         * Matrix.Band.t AdjointTypes'.DirectTypes.linsys_fn
  | BDlsBandCallbackSens
      of Matrix.Band.t AdjointTypes'.DirectTypes.jac_callback_with_sens
         * Matrix.Band.t AdjointTypes'.DirectTypes.linsys_fn

  (* Sls *)
  | SlsKluCallback
      : ('s Matrix.Sparse.t) DirectTypes.jac_callback
        * ('s Matrix.Sparse.t) DirectTypes.linsys_fn
        -> ('a, 'kind) linsolv_callbacks
  | BSlsKluCallback
      : ('s Matrix.Sparse.t) AdjointTypes'.DirectTypes.jac_callback_no_sens
         * ('s Matrix.Sparse.t) AdjointTypes'.DirectTypes.linsys_fn
        -> ('a, 'kind) linsolv_callbacks
  | BSlsKluCallbackSens
      : ('s Matrix.Sparse.t) AdjointTypes'.DirectTypes.jac_callback_with_sens
         * ('s Matrix.Sparse.t) AdjointTypes'.DirectTypes.linsys_fn
        -> ('a, 'kind) linsolv_callbacks

  | SlsSuperlumtCallback
      : ('s Matrix.Sparse.t) DirectTypes.jac_callback
        * ('s Matrix.Sparse.t) DirectTypes.linsys_fn
        -> ('a, 'kind) linsolv_callbacks
  | BSlsSuperlumtCallback
      : ('s Matrix.Sparse.t) AdjointTypes'.DirectTypes.jac_callback_no_sens
         * ('s Matrix.Sparse.t) AdjointTypes'.DirectTypes.linsys_fn
        -> ('a, 'kind) linsolv_callbacks
  | BSlsSuperlumtCallbackSens
      : ('s Matrix.Sparse.t) AdjointTypes'.DirectTypes.jac_callback_with_sens
         * ('s Matrix.Sparse.t) AdjointTypes'.DirectTypes.linsys_fn
        -> ('a, 'kind) linsolv_callbacks

  (* Custom *)
  | DirectCustomCallback
      : 'm DirectTypes.jac_callback
        * 'm DirectTypes.linsys_fn
        -> ('a, 'kind) linsolv_callbacks
  | BDirectCustomCallback
      : 'm AdjointTypes'.DirectTypes.jac_callback_no_sens
         * 'm AdjointTypes'.DirectTypes.linsys_fn
      -> ('a, 'kind) linsolv_callbacks
  | BDirectCustomCallbackSens
      : 'm AdjointTypes'.DirectTypes.jac_callback_with_sens
         * 'm AdjointTypes'.DirectTypes.linsys_fn
      -> ('a, 'kind) linsolv_callbacks

  (* Spils *)
  | SpilsCallback of 'a SpilsTypes'.jac_times_vec_fn option
                     * 'a SpilsTypes'.jac_times_setup_fn option
  | BSpilsCallback
      of 'a AdjointTypes'.SpilsTypes'.jac_times_vec_fn_no_sens option
         * 'a AdjointTypes'.SpilsTypes'.jac_times_setup_fn_no_sens option
  | BSpilsCallbackSens
      of 'a AdjointTypes'.SpilsTypes'.jac_times_vec_fn_with_sens option
         * 'a AdjointTypes'.SpilsTypes'.jac_times_setup_fn_with_sens option

and 'a linsolv_precfns =
  | NoPrecFns

  | PrecFns of 'a SpilsTypes'.precfns
  | BPrecFns of 'a AdjointTypes'.SpilsTypes'.precfns_no_sens
  | BPrecFnsSens of 'a AdjointTypes'.SpilsTypes'.precfns_with_sens

  | BandedPrecFns

  | BBDPrecFns of 'a CvodeBbdParamTypes.precfns
  | BBBDPrecFns of 'a CvodesBbdParamTypes.precfns

and ('a, 'kind) sensext =
    NoSensExt
  | FwdSensExt of ('a, 'kind) fsensext
  | BwdSensExt of ('a, 'kind) bsensext

and ('a, 'kind) fsensext = {
  (* Quadrature *)
  mutable quadrhsfn         : 'a QuadratureTypes.quadrhsfn;
  mutable checkquadvec      : (('a, 'kind) Nvector.t -> unit);
  mutable has_quad          : bool;

  (* Sensitivity *)
  mutable num_sensitivities : int;
  mutable sensarray1        : 'a array;
  mutable sensarray2        : 'a array;
  mutable senspvals         : RealArray.t option;
  (* keep a reference to prevent garbage collection *)

  mutable sensrhsfn         : 'a SensitivityTypes.sensrhsfn_all;
  mutable sensrhsfn1        : 'a SensitivityTypes.sensrhsfn1;
  mutable quadsensrhsfn     : 'a SensitivityTypes.QuadratureTypes.quadsensrhsfn;

  mutable fnls_solver       : ('a, 'kind, (('a, 'kind) session) NLSI.integrator)
                                NLSI.nonlinear_solver_hold;
  (* keep a reference to prevent garbage collection *)

  (* Adjoint *)
  mutable bsessions         : ('a, 'kind) session list;
  (* hold references to prevent garbage collection
     of backward sessions which are needed for
     callbacks. *)
}

and ('a, 'kind) bsensext = {
  (* Adjoint *)
  parent                : ('a, 'kind) session ;
  which                 : int;

  bnum_sensitivities    : int;
  bsensarray            : 'a array;

  mutable brhsfn          : 'a AdjointTypes'.brhsfn_no_sens;
  mutable brhsfn_sens     : 'a AdjointTypes'.brhsfn_with_sens;
  mutable bquadrhsfn      : 'a AdjointTypes'.QuadratureTypes.bquadrhsfn_no_sens;
  mutable bquadrhsfn_sens : 'a AdjointTypes'.QuadratureTypes.bquadrhsfn_with_sens;
  mutable checkbquadvec   : (('a, 'kind) Nvector.t -> unit);
}

(* Linear solver check functions *)

let ls_check_diag session =
  if Sundials_configuration.safe && session.ls_callbacks <> DiagNoCallbacks then
    raise LinearSolver.InvalidLinearSolver

let ls_check_direct session =
  if Sundials_configuration.safe then
    match session.ls_callbacks with
    | DlsDenseCallback _ | DlsBandCallback _
    | BDlsDenseCallback _ | BDlsDenseCallbackSens _
    | BDlsBandCallback _ | BDlsBandCallbackSens _
    | SlsKluCallback _ | BSlsKluCallback _ | BSlsKluCallbackSens _
    | SlsSuperlumtCallback _ | BSlsSuperlumtCallback _
    | BSlsSuperlumtCallbackSens _ | DirectCustomCallback _ -> ()
    | _ -> raise LinearSolver.InvalidLinearSolver

let ls_check_spils session =
  if Sundials_configuration.safe then
    match session.ls_callbacks with
    | SpilsCallback _ | BSpilsCallback _ | BSpilsCallbackSens _ -> ()
    | _ -> raise LinearSolver.InvalidLinearSolver

let ls_check_spils_band session =
  if Sundials_configuration.safe then
    match session.ls_precfns with
    | BandedPrecFns -> ()
    | _ -> raise LinearSolver.InvalidLinearSolver

let ls_check_spils_bbd session =
  if Sundials_configuration.safe then
    match session.ls_precfns with
    | BBDPrecFns _ | BBBDPrecFns _ -> ()
    | _ -> raise LinearSolver.InvalidLinearSolver

(* Types that depend on session *)

type 'kind serial_session = (Nvector_serial.data, 'kind) session
                            constraint 'kind = [>Nvector_serial.kind]

type ('data, 'kind) linear_solver =
  ('data, 'kind) session
  -> ('data, 'kind) Nvector.t
  -> unit

type 'k serial_linear_solver =
  (Nvector_serial.data, 'k) linear_solver
  constraint 'k = [>Nvector_serial.kind]

module SpilsTypes = struct
  include SpilsTypes'

  type ('a, 'k) set_preconditioner =
    ('a, 'k) session
    -> ('a, 'k) Nvector.t
    -> unit

  type ('a, 'k) preconditioner =
    LSI.Iterative.preconditioning_type * ('a, 'k) set_preconditioner

  type 'k serial_preconditioner = (Nvector_serial.data, 'k) preconditioner
                                  constraint 'k = [>Nvector_serial.kind]

end

module AdjointTypes = struct
  include AdjointTypes'
  (* Backwards session. *)
  type ('a, 'k) bsession = Bsession of ('a, 'k) session
  type 'k serial_bsession = (Nvector_serial.data, 'k) bsession
                            constraint 'k = [>Nvector_serial.kind]
  let tosession (Bsession s) = s
  let parent_and_which s =
    match (tosession s).sensext with
    | BwdSensExt se -> (se.parent, se.which)
    | _ -> failwith "Internal error: bsession invalid"

  type ('data, 'kind) linear_solver =
    ('data, 'kind) bsession
    -> ('data, 'kind) nvector
    -> unit
  type 'kind serial_linear_solver =
    (Nvector_serial.data, 'kind) linear_solver
    constraint 'kind = [>Nvector_serial.kind]

  module SpilsTypes = struct
    include SpilsTypes'

    type ('a, 'k) set_preconditioner =
      ('a, 'k) bsession -> ('a, 'k) session -> int -> ('a, 'k) nvector -> unit

    type ('a, 'k) preconditioner =
      LSI.Iterative.preconditioning_type * ('a, 'k) set_preconditioner

    type 'k serial_preconditioner = (Nvector_serial.data, 'k) preconditioner
                                    constraint 'k = [>Nvector_serial.kind]
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
let dummy_rhsfn _ _ _ =
  Sundials_impl.crash "Internal error: dummy_resfn called\n"
let dummy_rootsfn _ _ _ =
  Sundials_impl.crash "Internal error: dummy_rootsfn called\n"
let dummy_errh _ =
  Sundials_impl.crash "Internal error: dummy_errh called\n"
let dummy_errw _ _ =
  Sundials_impl.crash "Internal error: dummy_errw called\n"
let dummy_projfn _ _ _ _ _ =
  Sundials_impl.crash "Internal error: dummy_projfn called\n"
let dummy_monitorfn _ =
  Sundials_impl.crash "Internal error: dummy_monitorfn called\n"
let dummy_brhsfn_no_sens _ _ =
  Sundials_impl.crash "Internal error: dummy_brhsfn_no_sens called\n"
let dummy_brhsfn_with_sens _ _ _ =
  Sundials_impl.crash "Internal error: dummy_brhsfn_with_sens called\n"
let dummy_bquadrhsfn_no_sens _ _ =
  Sundials_impl.crash "Internal error: dummy_bquadrhsfn_no_sens called\n"
let dummy_bquadrhsfn_with_sens _ _ _ =
  Sundials_impl.crash "Internal error: dummy_bquadrhsfn_with_sens called\n"
let dummy_quadrhsfn _ _ _ =
  Sundials_impl.crash "Internal error: dummy_quadrhsfn called\n"
let dummy_sensrhsfn _ _ _ =
  Sundials_impl.crash "Internal error: dummy_sensresfn called\n"
let dummy_sensrhsfn1 _ _ _ _ =
  Sundials_impl.crash "Internal error: dummy_sensresfn1 called\n"
let dummy_quadsensrhsfn _ _ =
  Sundials_impl.crash "Internal error: dummy_quadsensrhsfn called\n"
