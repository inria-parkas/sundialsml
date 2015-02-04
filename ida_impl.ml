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

(* Types shared between Ida, Idas, Ida_bbd, and Idas_bbd.  See the
   notes on Cvode_impl about the rationale behind this module.  *)

(*
 * NB: The order of variant constructors and record fields is important!
 *     If these types are changed or augmented, the corresponding declarations
 *     in cvode_ml.h (and code in cvode_ml.c) must also be updated.
 *)

(* Dummy callbacks.  These dummes getting called indicates a fatal
   bug.  Rather than raise an exception (which may or may not get
   propagated properly depending on the context), we immediately abort
   the program. *)
external crash : string -> unit = "sundials_crash"

type 'a double = 'a * 'a
type 'a triple = 'a * 'a * 'a

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
    (Sundials.RealArray.t triple, Sundials.RealArray.t) jacobian_arg
    -> Dls.DenseMatrix.t
    -> unit

  (* These fields are accessed from cvode_ml.c *)
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
    -> (Sundials.RealArray.t triple, Sundials.RealArray.t) jacobian_arg
    -> Dls.BandMatrix.t
    -> unit

  (* These fields are accessed from cvode_ml.c *)
  type band_jac_callback =
    {
      bjacfn: band_jac_fn;
      mutable bmat : Dls.BandMatrix.t option
    }

  let no_band_callback = {
      bjacfn = (fun _ _ _ -> crash "no band callback");
      bmat = None;
    }
end

module SpilsCommonTypes = struct
  (* Types that don't depend on jacobian_arg.  *)
  type gramschmidt_type = Spils.gramschmidt_type =
    | ModifiedGS
    | ClassicalGS
end

module SpilsTypes' = struct
  include SpilsCommonTypes
  type 'a prec_solve_fn =
    ('a, 'a) jacobian_arg
    -> 'a
    -> 'a
    -> float
    -> unit
  type 'a prec_setup_fn = ('a triple, 'a) jacobian_arg -> unit
  type 'a jac_times_vec_fn =
    ('a double, 'a) jacobian_arg
    -> 'a           (* v *)
    -> 'a           (* Jv *)
    -> unit

  type 'a callbacks =
    {
      prec_solve_fn : 'a prec_solve_fn;
      prec_setup_fn : 'a prec_setup_fn option;
      jac_times_vec_fn : 'a jac_times_vec_fn option;
    }

  let no_prec_callbacks = {
      prec_solve_fn = (fun _ _ _ _ -> crash "no prec solve callback");
      prec_setup_fn = None;
      jac_times_vec_fn = None;
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
  type 'd sensresfn_args =
    {
      t : float;
      y : 'd;
      y' : 'd;
      res : 'd;
      s : 'd array;
      s' : 'd array;
      tmp : 'd triple
    }

  type 'd sensresfn = 'd sensresfn_args -> 'd array -> unit

  module QuadratureTypes = struct
    type 'd quadsensrhsfn =
      'd quadsensrhsfn_args
      -> 'd array
      -> unit

    and 'd quadsensrhsfn_args =
      {
        t : float;
        y : 'd;
        y' : 'd;
        s : 'd array;
        s' : 'd array;
        q : 'd;
        tmp : 'd triple;
      }
  end
end

module AdjointTypes' = struct
  type 'd bresfn_args =
    {
      t : float;
      y : 'd;
      y' : 'd;
      yb : 'd;
      yb' : 'd;
    }
  type 'a bresfn_no_sens = 'a bresfn_args -> 'a -> unit
  and 'a bresfn_with_sens = 'a bresfn_args -> 'a array -> 'a array -> 'a -> unit

  type 'a bresfn =
      NoSens of 'a bresfn_no_sens
    | WithSens of 'a bresfn_with_sens

  module QuadratureTypes = struct
    type 'd bquadrhsfn_args =
      {
        t : float;
        y : 'd;
        y' : 'd;
        yb : 'd;
        yb' : 'd;
      }

    type 'a bquadrhsfn =
        NoSens of 'a bquadrhsfn_no_sens
      | WithSens of 'a bquadrhsfn_with_sens

    and 'a bquadrhsfn_no_sens = 'a bquadrhsfn_args -> 'a -> unit
    and 'a bquadrhsfn_with_sens = 'a bquadrhsfn_args -> 'a array -> 'a array ->
                                  'a -> unit
  end

  type ('t, 'a) jacobian_arg =
    {
      jac_t   : float;
      jac_y   : 'a;
      jac_y'  : 'a;
      jac_yb  : 'a;
      jac_yb' : 'a;
      jac_resb : 'a;
      jac_coef : float;
      jac_tmp : 't
    }

  (* This is NOT the same as DlsTypes defined above.  This version
     refers to a different jacobian_arg, the one that was just
     defined.  *)

  module DlsTypes = struct
    type dense_jac_fn =
      (Sundials.RealArray.t triple, Sundials.RealArray.t) jacobian_arg
      -> Dls.DenseMatrix.t
      -> unit

    (* These fields are accessed from cvode_ml.c *)
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
      -> (Sundials.RealArray.t triple, Sundials.RealArray.t) jacobian_arg
      -> Dls.BandMatrix.t
      -> unit

    (* These fields are accessed from cvode_ml.c *)
    type band_jac_callback =
      {
        bjacfn: band_jac_fn;
        mutable bmat : Dls.BandMatrix.t option
      }

    let no_band_callback = {
        bjacfn = (fun _ _ _ -> crash "no band callback");
        bmat = None;
      }
  end

  (* Ditto. *)
  module SpilsTypes' = struct
    include SpilsCommonTypes
    type 'a prec_solve_fn =
      ('a, 'a) jacobian_arg
      -> 'a
      -> 'a
      -> float
      -> unit
    type 'a prec_setup_fn = ('a triple, 'a) jacobian_arg -> unit
    type 'a jac_times_vec_fn =
      ('a, 'a) jacobian_arg
      -> 'a
      -> 'a
      -> unit

    type 'a callbacks =
      {
        prec_solve_fn : 'a prec_solve_fn;
        prec_setup_fn : 'a prec_setup_fn option;
        jac_times_vec_fn : 'a jac_times_vec_fn option;
      }

    let no_prec_callbacks = {
        prec_solve_fn = (fun _ _ _ _ -> crash "no prec solve callback");
        prec_setup_fn = None;
        jac_times_vec_fn = None;
      }
  end
end

module IdasBbdParamTypes = struct
  type 'a local_fn = 'a AdjointTypes'.bresfn_args -> 'a -> unit
  type 'a comm_fn = 'a AdjointTypes'.bresfn_args -> unit
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
type 'a rootsfn = float -> 'a -> 'a -> Sundials.RealArray.t -> unit
type error_handler = Sundials.error_details -> unit
type 'a error_weight_fun = 'a -> 'a -> unit

(* Session: here comes the big blob.  These mutually recursive types
   cannot be handed out separately to modules without menial
   repetition, so we'll just have them all here, at the top of the
   Types module.  *)

type ('a,'kind) session = {
  ida        : ida_mem;
  backref    : c_weak_ref;
  nroots     : int;
  err_file   : ida_file;
  checkvec   : (('a, 'kind) Nvector.t -> unit);

  (* Temporary storage for exceptions raised within callbacks.  *)
  mutable exn_temp   : exn option;
  (* Tracks whether IDASetId has been called. *)
  mutable id_set     : bool;

  mutable resfn      : 'a resfn;
  mutable rootsfn    : 'a rootsfn;
  mutable errh       : error_handler;
  mutable errw       : 'a error_weight_fun;

  mutable ls_callbacks : ('a, 'kind) linsolv_callbacks;

  mutable sensext      : ('a, 'kind) sensext; (* Used by IDAS *)
}

and ('a, 'kind) sensext =
    NoSensExt
  | FwdSensExt of ('a, 'kind) fsensext
  | BwdSensExt of ('a, 'kind) bsensext

and ('a, 'kind) fsensext = {
  (* Quadrature *)
  mutable quadrhsfn       : 'a QuadratureTypes.quadrhsfn;
  mutable checkquadvec    : (('a, 'kind) Nvector.t -> unit);
  mutable has_quad        : bool;

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

  mutable bresfn          : 'a AdjointTypes'.bresfn_no_sens;
  mutable bresfn_sens     : 'a AdjointTypes'.bresfn_with_sens;
  mutable bquadrhsfn      : 'a AdjointTypes'.QuadratureTypes.bquadrhsfn_no_sens;
  mutable bquadrhsfn_sens : 'a AdjointTypes'.QuadratureTypes.bquadrhsfn_with_sens;
  mutable checkbquadvec   : (('a, 'kind) Nvector.t -> unit);
}

and ('a, 'kind) linsolv_callbacks =
  | NoCallbacks

  (* Dls *)
  | DlsDenseCallback of DlsTypes.dense_jac_callback
  | DlsBandCallback  of DlsTypes.band_jac_callback

  | BDlsDenseCallback of AdjointTypes'.DlsTypes.dense_jac_callback
  | BDlsBandCallback  of AdjointTypes'.DlsTypes.band_jac_callback

  (* Spils *)

  | SpilsCallback of 'a SpilsTypes'.callbacks
  | SpilsBBDCallback of 'a IdaBbdParamTypes.callbacks

  | BSpilsCallback of 'a AdjointTypes'.SpilsTypes'.callbacks
  | BSpilsBBDCallback of 'a IdasBbdParamTypes.callbacks

  (* Alternate *)
  | AlternateCallback of ('a, 'kind) alternate_linsolv

and ('data, 'kind) alternate_linsolv =
  {
    linit  : ('data, 'kind) linit' option;
    lsetup : ('data, 'kind) lsetup' option;
    lsolve : ('data, 'kind) lsolve';
  }
and ('data, 'kind) linit' = ('data, 'kind) session -> unit
and 'data alternate_lsetup_args =
  {
    lsetup_y : 'data;
    lsetup_y' : 'data;
    lsetup_res : 'data;
    lsetup_tmp : 'data triple;
  }
and ('data, 'kind) lsetup' =
  ('data, 'kind) session
  -> 'data alternate_lsetup_args
  -> unit
and 'data alternate_lsolve_args =
  {
    lsolve_ewt : 'data;
    lsolve_y : 'data;
    lsolve_y' : 'data;
    lsolve_res : 'data;
  }
and ('data, 'kind) lsolve' =
  ('data, 'kind) session
  -> 'data alternate_lsolve_args
  -> 'data
  -> unit

(* Linear solver check functions *)

let ls_check_dls session =
  if Sundials_config.safe then
    match session.ls_callbacks with
    | DlsDenseCallback _ | DlsBandCallback _
    | BDlsDenseCallback _ | BDlsBandCallback _ -> ()
    | _ -> raise Sundials.InvalidLinearSolver

let ls_check_spils session =
  if Sundials_config.safe then
    match session.ls_callbacks with
    | SpilsCallback _ | SpilsBBDCallback _
    | BSpilsCallback _ | BSpilsBBDCallback _ -> ()
    | _ -> raise Sundials.InvalidLinearSolver

let ls_check_spils_bbd session =
  if Sundials_config.safe then
    match session.ls_callbacks with
    | SpilsBBDCallback _ | BSpilsBBDCallback _ -> ()
    | _ -> raise Sundials.InvalidLinearSolver

(* Types that depend on session *)

type serial_session = (Nvector_serial.data, Nvector_serial.kind) session

(* IDA's linear_solver receives two vectors, y and y'.  They usually
   (always?) have identical size and other properties, so one of them
   can be safely ignored.  It's just in case that both are given.  *)
type ('data, 'kind) linear_solver =
  ('data, 'kind) session
  -> ('data, 'kind) Nvector.t (* y *)
  -> ('data, 'kind) Nvector.t (* y' *)
  -> unit

type serial_linear_solver =
  (Nvector_serial.data, Nvector_serial.kind) linear_solver

module SpilsTypes = struct
  include SpilsTypes'

  type ('a, 'k) set_preconditioner =
    ('a, 'k) session -> ('a, 'k) Nvector.t -> ('a, 'k) Nvector.t -> unit

  (* IDA(S) supports only left preconditioning.  *)
  type ('a, 'k) preconditioner =
    | InternalPrecNone of ('a, 'k) set_preconditioner
    | InternalPrecLeft of ('a, 'k) set_preconditioner

  type serial_preconditioner =
    (Nvector_serial.data, Nvector_serial.kind) preconditioner

end

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
  and 'data lsetup_args = 'data alternate_lsetup_args =
    {
      lsetup_y : 'data;
      lsetup_y' : 'data;
      lsetup_res : 'data;
      lsetup_tmp : 'data triple;
    }
  and 'data lsolve_args = 'data alternate_lsolve_args =
    {
      lsolve_ewt : 'data;
      lsolve_y : 'data;
      lsolve_y' : 'data;
      lsolve_res : 'data;
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
    -> ('data, 'kind) Nvector.t (* y *)
    -> ('data, 'kind) Nvector.t (* y' *)
    -> unit
  type serial_linear_solver =
    (Nvector_serial.data, Nvector_serial.kind) linear_solver

  module SpilsTypes = struct
    include SpilsTypes'

    type ('a, 'k) set_preconditioner =
      ('a, 'k) bsession -> ('a, 'k) session -> int ->
      ('a, 'k) Nvector.t -> ('a, 'k) Nvector.t -> unit

    (* IDA(S) supports only left preconditioning.  *)
    type ('a, 'k) preconditioner =
      | InternalPrecNone of ('a, 'k) set_preconditioner
      | InternalPrecLeft of ('a, 'k) set_preconditioner

    type serial_preconditioner =
      (Nvector_serial.data, Nvector_serial.kind) preconditioner

  end
end

let read_weak_ref x : ('a, 'kind) session =
  match Weak.get x 0 with
  | Some y -> y
  | None -> raise (Failure "Internal error: weak reference is dead")

let dummy_resfn _ _ _ _ =
  crash "Internal error: dummy_resfn called\n"
let dummy_rootsfn _ _ _ _ =
  crash "Internal error: dummy_rootsfn called\n"
let dummy_errh _ =
  crash "Internal error: dummy_errh called\n"
let dummy_errw _ _ =
  crash "Internal error: dummy_errw called\n"
let dummy_bresfn_no_sens _ _ =
  crash "Internal error: dummy_bresfn_no_sens called\n"
let dummy_bresfn_with_sens _ _ _ _ =
  crash "Internal error: dummy_bresfn_with_sens called\n"
let dummy_bquadrhsfn_no_sens _ _ =
  crash "Internal error: dummy_bquadrhsfn_no_sens called\n"
let dummy_bquadrhsfn_with_sens _ _ _ _ =
  crash "Internal error: dummy_bquadrhsfn_with_sens called\n"
let dummy_quadrhsfn _ _ _ _ =
  crash "Internal error: dummy_quadrhsfn called\n"
let dummy_sensresfn _ _ =
  crash "Internal error: dummy_sensresfn called\n"
let dummy_quadsensrhsfn _ _ =
  crash "Internal error: dummy_quadsensrhsfn called\n"
