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

(* Types shared between Cvode, Cvodes, Cvode_bbd, and Cvodes_bbd.  *)

(* This module's purpose is just to define types shared between Cvode,
   Cvodes, Cvode_bbd, and Cvode_bbds.  However, the session and
   linear_solver types must have their implementation exposed to all
   four modules yet be abstract to code that lives outside of the
   sundialsml library.

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

(*
 * NB: The order of variant constructors and record fields is important!
 *     If these types are changed or augmented, the corresponding declarations
 *     in cvode_ml.h (and code in cvode_ml.c) must also be updated.
 *)

type ('data, 'kind) nvector = ('data, 'kind) Nvector.t
module RealArray = Sundials.RealArray

type 'a double = 'a * 'a
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

  (* These fields are accessed from cvode_ml.c *)
  type dense_jac_callback =
    {
      jacfn: dense_jac_fn;
      mutable dmat : Dls.DenseMatrix.t option
    }

  type band_jac_fn =
    bandrange
    -> (RealArray.t triple, RealArray.t) jacobian_arg
    -> Dls.BandMatrix.t
    -> unit

  (* These fields are accessed from cvode_ml.c *)
  type band_jac_callback =
    {
      bjacfn: band_jac_fn;
      mutable bmat : Dls.BandMatrix.t option
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

  type 'a callbacks =
    {
      prec_solve_fn : 'a prec_solve_fn;
      prec_setup_fn : 'a prec_setup_fn option;
      jac_times_vec_fn : 'a jac_times_vec_fn option;
    }
end

module AlternateTypes' = struct
  type conv_fail =
    | NoFailures
    | FailBadJ
    | FailOther
end

module CvodeBbdParamTypes = struct
  type 'a local_fn = float -> 'a -> 'a -> unit
  type 'a comm_fn = float -> 'a -> unit
  type 'a callbacks =
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
        yS : 'd array;
        yQ' : 'd;
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
      yB : 'd;
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
        yB : 'd;
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
      jac_yB  : 'a;
      jac_fyB : 'a;
      jac_tmp : 't
    }

  (* This is NOT the same as DlsTypes defined above.  This version
     refers to a different jacobian_arg, the one that was just
     defined.  *)
  module DlsTypes = struct
    type dense_jac_fn =
      (RealArray.t triple, RealArray.t) jacobian_arg
      -> Dls.DenseMatrix.t
      -> unit

    (* These fields are accessed from cvodes_ml.c *)
    type dense_jac_callback =
      {
        jacfn: dense_jac_fn;
        mutable dmat : Dls.DenseMatrix.t option
      }

    type band_jac_fn =
      bandrange
      -> (RealArray.t triple, RealArray.t) jacobian_arg
      -> Dls.BandMatrix.t
      -> unit

    (* These fields are accessed from cvodes_ml.c *)
    type band_jac_callback =
      {
        bjacfn: band_jac_fn;
        mutable bmat : Dls.BandMatrix.t option
      }
  end

  (* Ditto. *)
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

    type 'a callbacks =
      {
        prec_solve_fn : 'a prec_solve_fn;
        prec_setup_fn : 'a prec_setup_fn option;
        jac_times_vec_fn : 'a jac_times_vec_fn option;
      }
  end
end

module CvodesBbdParamTypes = struct
  type 'a local_fn = float -> 'a -> 'a -> 'a -> unit
  type 'a comm_fn = float -> 'a -> 'a -> unit
  type 'a callbacks =
    {
      local_fn : 'a local_fn;
      comm_fn  : 'a comm_fn option;
    }
end
module CvodesBbdTypes = CvodeBbdTypes

type cvode_mem
type cvode_file
type c_weak_ref

type 'a rhsfn = float -> 'a -> 'a -> unit
type 'a rootsfn = float -> 'a -> Sundials.RealArray.t -> unit
type error_handler = Sundials.error_details -> unit
type 'a error_weight_fun = 'a -> 'a -> unit

(* Session: here comes the big blob.  These mutually recursive types
   cannot be handed out separately to modules without menial
   repetition, so we'll just have them all here, at the top of the
   Types module.  *)

type ('a, 'kind) session = {
  cvode      : cvode_mem;
  backref    : c_weak_ref;
  nroots     : int;
  err_file   : cvode_file;
  checkvec   : (('a, 'kind) Nvector.t -> unit);

  mutable exn_temp     : exn option;

  mutable rhsfn        : 'a rhsfn;
  mutable rootsfn      : 'a rootsfn;
  mutable errh         : error_handler;
  mutable errw         : 'a error_weight_fun;

  mutable ls_callbacks : ('a, 'kind) linsolv_callbacks;

  mutable sensext      : ('a, 'kind) sensext (* Used by Cvodes *)
}

and ('a, 'kind) linsolv_callbacks =
  | NoCallbacks

  | DenseCallback of DlsTypes.dense_jac_callback
  | BandCallback  of DlsTypes.band_jac_callback
  | SpilsCallback of 'a SpilsTypes'.callbacks
  | SpilsBandCallback of 'a SpilsTypes'.jac_times_vec_fn option
                         (* Invariant: 'a = RealArray.t *)
  | BBDCallback of 'a CvodeBbdParamTypes.callbacks

  | AlternateCallback of ('a, 'kind) alternate_linsolv

  | BDenseCallback of AdjointTypes'.DlsTypes.dense_jac_callback
  | BBandCallback  of AdjointTypes'.DlsTypes.band_jac_callback
  | BSpilsCallback of 'a AdjointTypes'.SpilsTypes'.callbacks
  | BSpilsBandCallback of 'a AdjointTypes'.SpilsTypes'.jac_times_vec_fn option
                          (* Invariant: 'a = RealArray.t *)
  | BBBDCallback of 'a CvodesBbdParamTypes.callbacks

and ('a, 'kind) sensext =
    NoSensExt
  | FwdSensExt of ('a, 'kind) fsensext
  | BwdSensExt of ('a, 'kind) bsensext

and ('a, 'kind) fsensext = {
  (* Quadrature *)
  mutable quadrhsfn         : 'a QuadratureTypes.quadrhsfn;
  mutable checkquadvec      : (('a, 'kind) Nvector.t -> unit);

  (* Sensitivity *)
  mutable num_sensitivities : int;
  mutable sensarray1        : 'a array;
  mutable sensarray2        : 'a array;
  mutable senspvals         : RealArray.t option;
  (* keep a reference to prevent garbage collection *)

  mutable sensrhsfn         : 'a SensitivityTypes.sensrhsfn_all;
  mutable sensrhsfn1        : 'a SensitivityTypes.sensrhsfn1;
  mutable quadsensrhsfn     : 'a SensitivityTypes.QuadratureTypes.quadsensrhsfn;
  mutable checkquadsensvec  : (('a, 'kind) Nvector.t -> unit);

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

(* Types that depend on session *)

type serial_session = (Nvector_serial.data, Nvector_serial.kind) session

type ('data, 'kind) linear_solver =
  ('data, 'kind) session
  -> ('data, 'kind) nvector
  -> unit

type serial_linear_solver =
  (Nvector_serial.data, Nvector_serial.kind) linear_solver

module SpilsTypes = struct
  include SpilsTypes'

  type ('a, 'k) set_preconditioner =
    ('a, 'k) session -> ('a, 'k) nvector -> unit

  type ('a, 'k) preconditioner =
    | InternalPrecNone
    | InternalPrecLeft of ('a, 'k) set_preconditioner
    | InternalPrecRight of ('a, 'k) set_preconditioner
    | InternalPrecBoth of ('a, 'k) set_preconditioner

  type serial_preconditioner =
    (Nvector_serial.data, Nvector_serial.kind) preconditioner
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

end

module AdjointTypes = struct
  include AdjointTypes'
  (* Backwards session. *)
  type ('a, 'k) bsession = Bsession of ('a, 'k) session
  type serial_bsession = (Nvector_serial.data, Nvector_serial.kind) bsession
  let tosession (Bsession s) = s

  type ('data, 'kind) linear_solver =
    ('data, 'kind) bsession
    -> ('data, 'kind) nvector
    -> unit
  type serial_linear_solver =
    (Nvector_serial.data, Nvector_serial.kind) linear_solver

  module SpilsTypes = struct
    include SpilsTypes'

    type ('a, 'k) set_preconditioner =
      ('a, 'k) bsession -> ('a, 'k) session -> int -> ('a, 'k) nvector -> unit

    type ('a, 'k) preconditioner =
      | InternalPrecNone
      | InternalPrecLeft of ('a, 'k) set_preconditioner
      | InternalPrecRight of ('a, 'k) set_preconditioner
      | InternalPrecBoth of ('a, 'k) set_preconditioner

    type serial_preconditioner =
      (Nvector_serial.data, Nvector_serial.kind) preconditioner
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
external crash : string -> unit = "sundials_crash"
let dummy_rhsfn _ _ _ =
  crash "Internal error: dummy_resfn called\n"
let dummy_rootsfn _ _ _ =
  crash "Internal error: dummy_rootsfn called\n"
let dummy_errh _ =
  crash "Internal error: dummy_errh called\n"
let dummy_errw _ _ =
  crash "Internal error: dummy_errw called\n"
let dummy_brhsfn_no_sens _ _ =
  crash "Internal error: dummy_brhsfn_no_sens called\n"
let dummy_brhsfn_with_sens _ _ _ =
  crash "Internal error: dummy_brhsfn_with_sens called\n"
let dummy_bquadrhsfn_no_sens _ _ =
  crash "Internal error: dummy_bquadrhsfn_no_sens called\n"
let dummy_bquadrhsfn_with_sens _ _ _ =
  crash "Internal error: dummy_bquadrhsfn_with_sens called\n"
let dummy_quadrhsfn _ _ _ =
  crash "Internal error: dummy_quadrhsfn called\n"
let dummy_sensrhsfn _ _ _ =
  crash "Internal error: dummy_sensresfn called\n"
let dummy_sensrhsfn1 _ _ _ _ =
  crash "Internal error: dummy_sensresfn called\n"
let dummy_quadsensrhsfn _ _ =
  crash "Internal error: dummy_quadsensrhsfn called\n"
