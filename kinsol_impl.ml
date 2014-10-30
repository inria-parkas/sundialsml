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

(* Types shared between Kinsol and Kinsol_bbd.  See the notes on
   Cvode_impl about the rationale behind this module.  *)

(*
 * NB: The order of variant constructors and record fields is important!
 *     If these types are changed or augmented, the corresponding declarations
 *     in cvode_ml.h (and code in cvode_ml.c) must also be updated.
 *)

type ('a, 'k) nvector = ('a, 'k) Nvector.t

type 'a single_tmp = 'a
type 'a double_tmp = 'a * 'a

type ('t, 'a) jacobian_arg =
  {
    jac_u   : 'a;
    jac_fu  : 'a;
    jac_tmp : 't
  }

type real_array = Sundials.RealArray.t

type bandrange = { mupper : int; mlower : int; }

module DlsTypes = struct
  type dense_jac_fn =
    (real_array double_tmp, real_array) jacobian_arg
    -> Dls.DenseMatrix.t -> unit

  (* These fields are accessed from cvode_ml.c *)
  type dense_jac_callback =
    {
      jacfn: dense_jac_fn;
      mutable dmat : Dls.DenseMatrix.t option
    }

  type band_jac_fn =
    bandrange
    -> (real_array double_tmp, real_array) jacobian_arg
    -> Dls.BandMatrix.t
    -> unit

  (* These fields are accessed from cvode_ml.c *)
  type band_jac_callback =
    {
      bjacfn: band_jac_fn;
      mutable bmat : Dls.BandMatrix.t option
    }
end

module SpilsTypes' = struct
  type 'a solve_arg = { uscale : 'a; fscale : 'a; }

  type 'a prec_solve_fn =
    ('a single_tmp, 'a) jacobian_arg
    -> 'a solve_arg
    -> 'a
    -> unit

  type 'a prec_setup_fn =
    ('a double_tmp, 'a) jacobian_arg
    -> 'a solve_arg
    -> unit

  type 'a jac_times_vec_fn = 'a -> 'a -> 'a -> bool -> bool

  type 'a callbacks =
    {
      prec_solve_fn : 'a prec_solve_fn option;

      prec_setup_fn : 'a prec_setup_fn option;

      jac_times_vec_fn : 'a jac_times_vec_fn option;
    }
end

module KinsolBbdParamTypes = struct
  type 'a local_fn = 'a -> 'a -> unit
  type 'a comm_fn = 'a -> unit
  type 'a callbacks =
    {
      local_fn : 'a local_fn;
      comm_fn  : 'a comm_fn option;
    }
end

module KinsolBbdTypes = struct
  type bandwidths =
    {
      mudq    : int;
      mldq    : int;
      mukeep  : int;
      mlkeep  : int;
    }
end

type kin_mem
type kin_file
type c_weak_ref

type 'a sysfn = 'a -> 'a -> unit
type errh = Sundials.error_details -> unit
type infoh = Sundials.error_details -> unit

(* Session: here comes the big blob.  These mutually recursive types
   cannot be handed out separately to modules without menial
   repetition, so we'll just have them all here, at the top of the
   Types module.  *)

type ('a, 'k) session = {
  kinsol    : kin_mem;
  backref   : c_weak_ref;
  err_file  : kin_file;
  info_file : kin_file;
  checkvec  : (('a, 'k) Nvector.t -> unit);

  mutable neqs       : int;    (* only valid for 'kind = serial *)
  mutable exn_temp   : exn option;

  mutable sysfn      : 'a sysfn;
  mutable errh       : errh;
  mutable infoh      : infoh;

  mutable ls_callbacks : ('a, 'k) linsolv_callbacks;
}

and ('a, 'kind) linsolv_callbacks =
  | NoCallbacks

  | DenseCallback of DlsTypes.dense_jac_callback
  | BandCallback  of DlsTypes.band_jac_callback
  | SpilsCallback of 'a SpilsTypes'.callbacks

  | AlternateCallback of ('a, 'kind) alternate_linsolv

  | BBDCallback of 'a KinsolBbdParamTypes.callbacks

and ('data, 'kind) alternate_linsolv =
  {
    linit  : ('data, 'kind) linit' option;
    lsetup : ('data, 'kind) lsetup' option;
    lsolve : ('data, 'kind) lsolve';
  }
and ('data, 'kind) linit' = ('data, 'kind) session -> unit
and ('data, 'kind) lsetup' = ('data, 'kind) session -> unit
and ('data, 'kind) lsolve' =
  ('data, 'kind) session
  -> 'data
  -> 'data
  -> float option

(* Types that depend on session *)
type serial_session = (Nvector_serial.data, Nvector_serial.kind) session

type ('data, 'kind) linear_solver =
  ('data, 'kind) session
  -> ('data, 'kind) nvector option
  -> unit

module SpilsTypes = struct
  include SpilsTypes'

  type ('a, 'k) set_preconditioner =
    ('a, 'k) session -> ('a, 'k) nvector option -> unit
  
  type ('a, 'k) preconditioner =
    | InternalPrecNone
    | InternalPrecRight of ('a, 'k) set_preconditioner

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
end

let read_weak_ref x : ('a, 'k) session =
  match Weak.get x 0 with
  | Some y -> y
  | None -> raise (Failure "Internal error: weak reference is dead")

(* Dummy callbacks.  These dummes getting called indicates a fatal
   bug.  Rather than raise an exception (which may or may not get
   propagated properly depending on the context), we immediately abort
   the program. *)
external crash : string -> unit = "sundials_crash"

let dummy_sysfn _ _ =
  crash "Internal error: dummy_sysfn called\n"
let dummy_errh _ =
  crash "Internal error: dummy_errh called\n"
let dummy_infoh _ =
  crash "Internal error: dummy_infoh called\n"
