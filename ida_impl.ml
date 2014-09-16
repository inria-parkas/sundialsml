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

(* Types shared between Ida, Idas, Ida_bbd, and Idas_bbd.  See the
   notes on Cvode_impl about the rationale behind this module.  *)

(*
 * NB: The order of variant constructors and record fields is important!
 *     If these types are changed or augmented, the corresponding declarations
 *     in cvode_ml.h (and code in cvode_ml.c) must also be updated.
 *)

type ('data, 'kind) nvector = ('data, 'kind) Sundials.nvector

type 'a single_tmp = 'a
type 'a double_tmp = 'a * 'a
type 'a triple_tmp = 'a * 'a * 'a

type ('t, 'a) jacobian_arg =
  {
    jac_t    : float;
    jac_y    : 'a;
    jac_y'   : 'a;
    jac_res  : 'a;
    jac_coef : float;
    jac_tmp  : 't
  }

type bandrange = { mupper : int; mlower : int; }

module DlsTypes = struct
  type dense_jac_fn =
    (Sundials.RealArray.t triple_tmp, Sundials.RealArray.t) jacobian_arg
    -> Dls.DenseMatrix.t
    -> unit

  type band_jac_fn =
    bandrange
    -> (Sundials.RealArray.t triple_tmp, Sundials.RealArray.t) jacobian_arg
    -> Dls.BandMatrix.t
    -> unit
end

module SpilsCommonTypes = struct
  (* Types that don't depend on jacobian_arg.  *)
  type gramschmidt_type = Spils.gramschmidt_type =
    | ModifiedGS
    | ClassicalGS

  type preconditioning_type = Spils.preconditioning_type =
    | PrecNone
    | PrecLeft
    | PrecRight
    | PrecBoth
end

module SpilsTypes = struct
  include SpilsCommonTypes
  type 'a prec_solve_fn =
    ('a single_tmp, 'a) jacobian_arg
    -> 'a
    -> 'a
    -> float
    -> unit
  type 'a prec_setup_fn = ('a triple_tmp, 'a) jacobian_arg -> unit
  type 'a jac_times_vec_fn =
    ('a double_tmp, 'a) jacobian_arg
    -> 'a           (* v *)
    -> 'a           (* Jv *)
    -> unit

  type 'a callbacks =
    {
      prec_solve_fn : 'a prec_solve_fn option;
      prec_setup_fn : 'a prec_setup_fn option;
      jac_times_vec_fn : 'a jac_times_vec_fn option;
    }
end

module IdaBbdParamTypes = struct
  type 'a local_fn = float -> 'a -> 'a -> 'a  -> unit
  type 'a comm_fn = float -> 'a -> 'a -> unit
  type 'a callbacks =
    {
      local_fn : 'a local_fn;
      comm_fn  : 'a comm_fn option;
    }
end

module IdaBbdTypes = struct
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
  type 'a quadrhsfn = float -> 'a -> 'a -> 'a -> unit
end

module SensitivityTypes = struct
  type 'a sensresfn =
    float                            (* t *)
    -> 'a                            (* y *)
    -> 'a                            (* y' *)
    -> 'a                            (* resval *)
    -> 'a array                      (* yS *)
    -> 'a array                      (* y'S *)
    -> 'a array                      (* resvalS *)
    -> 'a                            (* tmp1 *)
    -> 'a                            (* tmp2 *)
    -> 'a                            (* tmp3 *)
    -> unit

  module QuadratureTypes = struct
    type 'a quadsensrhsfn =
      float          (* t *)
      -> 'a          (* y *)
      -> 'a          (* y' *)
      -> 'a array    (* yS *)
      -> 'a array    (* y'S *)
      -> 'a          (* rrQ *)
      -> 'a array    (* rhsvalQs *)
      -> 'a          (* tmp1 *)
      -> 'a          (* tmp2 *)
      -> 'a          (* tmp3 *)
      -> unit
  end
end

module AdjointTypes' = struct
  type 'a bresfn_no_sens =
    float             (* t *)
    -> 'a             (* y *)
    -> 'a             (* y' *)
    -> 'a             (* yB *)
    -> 'a             (* y'B *)
    -> 'a             (* resvalB *)
    -> unit
  and 'a bresfn_with_sens =
    float          (* t *)
    -> 'a          (* y *)
    -> 'a          (* y' *)
    -> 'a array    (* yS *)
    -> 'a array    (* y'S *)
    -> 'a          (* yB *)
    -> 'a          (* y'B *)
    -> 'a          (* resvalB *)
    -> unit
  type 'a bresfn =
      NoSens of 'a bresfn_no_sens
    | WithSens of 'a bresfn_with_sens

  module QuadratureTypes = struct
    type 'a bquadrhsfn =
        NoSens of 'a bquadrhsfn_no_sens
      | WithSens of 'a bquadrhsfn_with_sens

    and 'a bquadrhsfn_no_sens =
      float             (* t *)
      -> 'a             (* y *)
      -> 'a             (* y' *)
      -> 'a             (* yB *)
      -> 'a             (* y'B *)
      -> 'a             (* rhsvalBQS *)
      -> unit
    and 'a bquadrhsfn_with_sens =
      float          (* t *)
      -> 'a          (* y *)
      -> 'a          (* y' *)
      -> 'a array    (* yS *)
      -> 'a array    (* y'S *)
      -> 'a          (* yB *)
      -> 'a          (* y'B *)
      -> 'a          (* rhsvalBQS *)
      -> unit
  end

  type ('t, 'a) jacobian_arg =
    {
      jac_t   : float;
      jac_y   : 'a;
      jac_y'  : 'a;
      jac_yb  : 'a;
      jac_y'b : 'a;
      jac_resb : 'a;
      jac_coef : float;
      jac_tmp : 't
    }

  (* This is NOT the same as DlsTypes defined above.  This version
     refers to a different jacobian_arg, the one that was just
     defined.  *)

  module DlsTypes = struct
    type dense_jac_fn =
      (Sundials.RealArray.t triple_tmp, Sundials.RealArray.t) jacobian_arg
      -> Dls.DenseMatrix.t
      -> unit

    type band_jac_fn =
      bandrange
      -> (Sundials.RealArray.t triple_tmp, Sundials.RealArray.t) jacobian_arg
      -> Dls.BandMatrix.t
      -> unit
  end

  (* Ditto. *)
  module SpilsTypes = struct
    include SpilsCommonTypes
    type 'a prec_solve_fn =
      ('a single_tmp, 'a) jacobian_arg
      -> 'a
      -> 'a
      -> float
      -> unit
    type 'a prec_setup_fn = ('a triple_tmp, 'a) jacobian_arg -> unit
    type 'a jac_times_vec_fn =
      ('a single_tmp, 'a) jacobian_arg
      -> 'a
      -> 'a
      -> unit

    type 'a callbacks =
      {
        prec_solve_fn : 'a prec_solve_fn option;
        prec_setup_fn : 'a prec_setup_fn option;
        jac_times_vec_fn : 'a jac_times_vec_fn option;
      }
  end
end

module IdasBbdParamTypes = struct
  type 'a local_fn = float -> 'a -> 'a -> 'a -> 'a -> 'a -> unit
  type 'a comm_fn = float -> 'a -> 'a -> 'a -> 'a -> unit
  type 'a callbacks =
    {
      local_fn : 'a local_fn;
      comm_fn  : 'a comm_fn option;
    }
end
module IdasBbdTypes = IdaBbdTypes

type ida_mem
type c_weak_ref
type ida_file

type 'a resfn = float -> 'a -> 'a -> 'a -> unit
type 'a rootsfn = float -> 'a -> 'a -> Sundials.Roots.val_array -> unit
type errh = Sundials.error_details -> unit
type 'a errw = 'a -> 'a -> unit

(* Session: here comes the big blob.  These mutually recursive types
   cannot be handed out separately to modules without menial
   repetition, so we'll just have them all here, at the top of the
   Types module.  *)

type ('a,'kind) session = {
  ida        : ida_mem;
  backref    : c_weak_ref;
  nroots     : int;
  err_file   : ida_file;

  (* Temporary storage for exceptions raised within callbacks.  *)
  mutable exn_temp   : exn option;

  mutable resfn      : 'a resfn;
  mutable rootsfn    : 'a rootsfn;
  mutable errh       : Sundials.error_details -> unit;
  mutable errw       : 'a errw;

  mutable ls_callbacks : ('a, 'kind) linsolv_callbacks;

  mutable sensext      : ('a, 'kind) sensext; (* Used by IDAS *)

  (* To be manipulated from the C side only.  *)
  mutable safety_check_flags : int;
}
and ('a, 'kind) sensext =
    NoSensExt
  | FwdSensExt of ('a, 'kind) fsensext
  | BwdSensExt of ('a, 'kind) bsensext

and ('a, 'kind) fsensext = {
  (* Quadrature *)
  mutable quadrhsfn       : 'a QuadratureTypes.quadrhsfn;

  (* Sensitivity *)
  mutable num_sensitivities : int;
  mutable sensarray1        : 'a array;
  mutable sensarray2        : 'a array;
  mutable sensarray3        : 'a array;
  mutable senspvals         : Sundials.RealArray.t option;
  (* keep a reference to prevent garbage collection *)

  mutable sensresfn         : 'a SensitivityTypes.sensresfn;

  mutable quadsensrhsfn     : 'a SensitivityTypes.QuadratureTypes.quadsensrhsfn;

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

  (* FIXME: extract num_sensitivities from bsensarray1 *)
  bnum_sensitivities    : int;
  bsensarray1           : 'a array;
  bsensarray2           : 'a array;

  mutable resfnb        : 'a AdjointTypes'.bresfn_no_sens;
  mutable resfnbs       : 'a AdjointTypes'.bresfn_with_sens;
  mutable bquadrhsfn    : 'a AdjointTypes'.QuadratureTypes.bquadrhsfn_no_sens;
  mutable bquadrhsfn1   : 'a AdjointTypes'.QuadratureTypes.bquadrhsfn_with_sens;
}

and ('a, 'kind) linsolv_callbacks =
  | NoCallbacks

  | DenseCallback of DlsTypes.dense_jac_fn
  | BandCallback  of DlsTypes.band_jac_fn
  | SpilsCallback of 'a SpilsTypes.callbacks
  | BBDCallback of 'a IdaBbdParamTypes.callbacks

  | AlternateCallback of ('a, 'kind) alternate_linsolv

  | BDenseCallback of AdjointTypes'.DlsTypes.dense_jac_fn
  | BBandCallback  of AdjointTypes'.DlsTypes.band_jac_fn
  | BSpilsCallback of 'a AdjointTypes'.SpilsTypes.callbacks
  | BBBDCallback of 'a IdasBbdParamTypes.callbacks

and ('data, 'kind) alternate_linsolv =
  {
    linit  : ('data, 'kind) linit' option;
    lsetup : ('data, 'kind) lsetup' option;
    lsolve : ('data, 'kind) lsolve';
  }
and ('data, 'kind) linit' = ('data, 'kind) session -> unit
and ('data, 'kind) lsetup' =
  ('data, 'kind) session
  -> 'data
  -> 'data
  -> 'data
  -> 'data triple_tmp
  -> unit
and ('data, 'kind) lsolve' =
  ('data, 'kind) session
  -> 'data
  -> 'data
  -> 'data
  -> 'data
  -> 'data
  -> unit

(* Types that depend on session *)

type serial_session = (Nvector_serial.data, Nvector_serial.kind) session

(* IDA's linear_solver receives two vectors, y and y'.  They usually
   (always?) have identical size and other properties, so one of them
   can be safely ignored.  It's just in case that both are given.  *)
type ('data, 'kind) linear_solver =
  ('data, 'kind) session
  -> ('data, 'kind) nvector (* y *)
  -> ('data, 'kind) nvector (* y' *)
  -> unit

type serial_linear_solver =
  (Nvector_serial.data, Nvector_serial.kind) linear_solver

module AlternateTypes = struct
  type ('data, 'kind) callbacks = ('data, 'kind) alternate_linsolv =
    {
      linit  : ('data, 'kind) linit option;
      lsetup : ('data, 'kind) lsetup option;
      lsolve : ('data, 'kind) lsolve;
    }
  and ('data, 'kind) linit  = ('data, 'kind) linit'
  and ('data, 'kind) lsetup = ('data, 'kind) lsetup'
  and ('data, 'kind) lsolve = ('data, 'kind) lsolve'
end

module AdjointTypes = struct
  include AdjointTypes'
  (* Backwards session. *)
  type ('a, 'k) bsession = Bsession of ('a, 'k) session
  type serial_bsession = (Nvector_serial.data, Nvector_serial.kind) bsession
  let tosession (Bsession s) = s

  type ('data, 'kind) linear_solver =
    ('data, 'kind) bsession
    -> ('data, 'kind) nvector (* y *)
    -> ('data, 'kind) nvector (* y' *)
    -> unit
  type serial_linear_solver =
    (Nvector_serial.data, Nvector_serial.kind) linear_solver
end

let read_weak_ref x : ('a, 'kind) session =
  match Weak.get x 0 with
  | Some y -> y
  | None -> raise (Failure "Internal error: weak reference is dead")

let adjust_retcode = fun session check_recoverable f x ->
  try f x; 0
  with
  | Sundials.RecoverableFailure when check_recoverable -> 1
  | e -> (session.exn_temp <- Some e; -1)

(* Dummy callbacks.  These dummes getting called indicates a fatal
   bug.  Rather than raise an exception (which may or may not get
   propagated properly depending on the context), we immediately abort
   the program. *)
external crash : string -> unit = "sundials_crash"
let dummy_resfn _ _ _ _ =
  crash "Internal error: dummy_resfn called\n"
let dummy_rootsfn _ _ _ _ =
  crash "Internal error: dummy_rootsfn called\n"
let dummy_errh _ =
  crash "Internal error: dummy_errh called\n"
let dummy_errw _ _ =
  crash "Internal error: dummy_errw called\n"
let dummy_resfnb _ _ _ _ _ _ =
  crash "Internal error: dummy_resfnb called\n"
let dummy_resfnbs _ _ _ _ _ _ _ _ =
  crash "Internal error: dummy_resfnbs called\n"
let dummy_bquadrhsfn _ _ _ _ _ _ =
  crash "Internal error: dummy_bquadrhsfn called\n"
let dummy_bquadrhsfn1 _ _ _ _ _ _ _ _ =
  crash "Internal error: dummy_bquadrhsfn1 called\n"
let dummy_quadrhsfn _ _ _ _ =
  crash "Internal error: dummy_quadrhsfn called\n"
let dummy_sensresfn _ _ _ _ _ _ _ _ _ _ =
  crash "Internal error: dummy_sensresfn called\n"
let dummy_quadsensrhsfn _ _ _ _ _ _ _ _ _ _ =
  crash "Internal error: dummy_quadsensrhsfn called\n"
