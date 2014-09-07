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

(* basic ida types *)

type ('data, 'kind) nvector = ('data, 'kind) Sundials.nvector
type real_array = Sundials.RealArray.t

type root_array = Sundials.Roots.t
type root_val_array = Sundials.Roots.val_array

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

and 'a spils_callbacks =
  {
    prec_solve_fn : (('a single_tmp, 'a) jacobian_arg -> 'a -> 'a -> float
                     -> unit) option;
    prec_setup_fn : (('a triple_tmp, 'a) jacobian_arg -> unit) option;
    jac_times_vec_fn : (('a double_tmp, 'a) jacobian_arg
                        -> 'a           (* v *)
                        -> 'a           (* Jv *)
                        -> unit) option;
  }

type bandrange = { mupper : int; mlower : int; }

type dense_jac_fn = (real_array triple_tmp, real_array) jacobian_arg
                  -> Dls.DenseMatrix.t -> unit

type band_jac_fn = bandrange -> (real_array triple_tmp, real_array) jacobian_arg
                 -> Dls.BandMatrix.t -> unit


type 'a quadrhsfn = float -> 'a -> 'a -> 'a -> unit

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

(* BBD definitions *)
module Bbd =
  struct
    type 'data callbacks =
      {
        local_fn : float -> 'data -> 'data -> 'data  -> unit;
        comm_fn  : (float -> 'data -> 'data -> unit) option;
      }
  end

module B =
  struct
    type 'a resfnb = float             (* t *)
                     -> 'a             (* y *)
                     -> 'a             (* y' *)
                     -> 'a             (* yB *)
                     -> 'a             (* y'B *)
                     -> 'a             (* resvalB *)
                     -> unit
    and 'a resfnbs = float          (* t *)
                     -> 'a          (* y *)
                     -> 'a          (* y' *)
                     -> 'a array    (* yS *)
                     -> 'a array    (* y'S *)
                     -> 'a          (* yB *)
                     -> 'a          (* y'B *)
                     -> 'a          (* resvalB *)
                     -> unit
    type 'a bresfn =
            Basic of 'a resfnb
          | WithSens of 'a resfnbs

    type 'a bquadrhsfn =
            Basic of 'a bquadrhsfn_basic
          | WithSens of 'a bquadrhsfn_withsens
      (* FIXME: less ad-hoc names *)
    and 'a bquadrhsfn_basic =
      float             (* t *)
      -> 'a             (* y *)
      -> 'a             (* y' *)
      -> 'a             (* yB *)
      -> 'a             (* y'B *)
      -> 'a             (* rhsvalBQS *)
      -> unit
    and 'a bquadrhsfn_withsens =
      float          (* t *)
      -> 'a          (* y *)
      -> 'a          (* y' *)
      -> 'a array    (* yS *)
      -> 'a array    (* y'S *)
      -> 'a          (* yB *)
      -> 'a          (* y'B *)
      -> 'a          (* rhsvalBQS *)
      -> unit

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

    type 'a spils_callbacks =
      {
        prec_solve_fn : (('a single_tmp, 'a) jacobian_arg
                         -> 'a -> 'a -> float -> unit) option;
        prec_setup_fn : (('a triple_tmp, 'a) jacobian_arg -> unit) option;
        jac_times_vec_fn :
          (('a single_tmp, 'a) jacobian_arg
           -> 'a (* v *)
           -> 'a (* Jv *)
           -> unit) option;
      }

    type dense_jac_fn =
          (real_array triple_tmp, real_array) jacobian_arg
              -> Dls.DenseMatrix.t -> unit

    type band_jac_fn =
          bandrange -> (real_array triple_tmp, real_array) jacobian_arg
              -> Dls.BandMatrix.t -> unit

    module Bbd =
      struct
        type 'data callbacks =
          {
            local_fn : float -> 'data -> 'data -> 'data
                       -> 'data -> 'data -> unit;
            comm_fn  : (float -> 'data -> 'data -> 'data -> 'data -> unit)
                       option;
          }
      end
  end

(* the session type *)

type ida_mem
type c_weak_ref
type ida_file

type ('a,'kind) session = {
        ida        : ida_mem;
        backref    : c_weak_ref;
        nroots     : int;
        err_file   : ida_file;

        (* Temporary storage for exceptions raised within callbacks.  *)
        mutable exn_temp   : exn option;

        mutable resfn      : float -> 'a -> 'a -> 'a -> unit;
        mutable rootsfn    : float -> 'a -> 'a -> root_val_array -> unit;
        mutable errh       : Sundials.error_details -> unit;
        mutable errw       : 'a -> 'a -> unit;

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
    mutable quadrhsfn       : 'a quadrhsfn;

    (* Sensitivity *)
    mutable num_sensitivities : int;
    mutable sensarray1        : 'a array;
    mutable sensarray2        : 'a array;
    mutable sensarray3        : 'a array;
    mutable senspvals         : Sundials.RealArray.t option;
                            (* keep a reference to prevent garbage collection *)

    mutable sensresfn         : 'a sensresfn;

    mutable quadsensrhsfn     : 'a quadsensrhsfn;

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

    mutable resfnb        : 'a B.resfnb;
    mutable resfnbs       : 'a B.resfnbs;
    mutable bquadrhsfn    : 'a B.bquadrhsfn_basic;
    mutable bquadrhsfn1   : 'a B.bquadrhsfn_withsens;
  }

and ('a, 'kind) linsolv_callbacks =
  | NoCallbacks

  | DenseCallback of dense_jac_fn
  | BandCallback  of band_jac_fn
  | SpilsCallback of 'a spils_callbacks
  | BBDCallback of 'a Bbd.callbacks

  | AlternateCallback of ('a, 'kind) alternate_linsolv

  | BDenseCallback of B.dense_jac_fn
  | BBandCallback  of B.band_jac_fn
  | BSpilsCallback of 'a B.spils_callbacks
  | BBBDCallback of 'a B.Bbd.callbacks

and ('data, 'kind) alternate_linsolv =
  {
    linit  : (('data, 'kind) session -> unit) option;
    lsetup : (('data, 'kind) session -> 'data
              -> 'data -> 'data -> 'data triple_tmp -> unit) option;
    lsolve : ('data, 'kind) session -> 'data -> 'data
              -> 'data -> 'data -> 'data -> unit;
  }

type ('a, 'k) bsession = Bsession of ('a, 'k) session
let tosession = function Bsession s -> s

(* IDA's linear_solver receives two vectors, y and y'.  They usually
   (always?) have identical size and other properties, so one of them
   can be safely ignored.  *)
type ('data, 'kind) linear_solver = ('data, 'kind) session
                                    -> ('data, 'kind) nvector (* y *)
                                    -> ('data, 'kind) nvector (* y' *)
                                    -> unit
type ('data, 'kind) blinear_solver = ('data, 'kind) bsession
                                     -> ('data, 'kind) nvector (* y *)
                                     -> ('data, 'kind) nvector (* y' *)
                                     -> unit

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
