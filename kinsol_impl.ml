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

type ('a, 'k) nvector = ('a, 'k) Sundials.nvector

type 'a single_tmp = 'a
type 'a double_tmp = 'a * 'a

type ('t, 'a) jacobian_arg =
  {
    jac_u   : 'a;
    jac_fu  : 'a;
    jac_tmp : 't
  }

type real_array = Sundials.RealArray.t

type dense_jac_fn = (real_array double_tmp, real_array) jacobian_arg
                        -> Dls.DenseMatrix.t -> unit

type bandrange = { mupper : int; mlower : int; }

type band_jac_fn = bandrange -> (real_array double_tmp, real_array) jacobian_arg
                        -> Dls.BandMatrix.t -> unit

type 'a prec_solve_arg = { uscale : 'a; fscale : 'a; }

type 'a spils_callbacks =
  {
    prec_solve_fn : (('a single_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                      -> 'a -> unit) option;

    prec_setup_fn : (('a double_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                     -> unit) option;

    jac_times_vec_fn : ('a -> 'a -> 'a -> bool -> bool) option;
  }

(* BBD definitions *)
module Bbd =
  struct
    type 'data callbacks =
      {
        local_fn : 'data -> 'data -> unit;
        comm_fn  : ('data -> unit) option;
      }
  end

type kin_mem
type kin_file
type c_weak_ref

type ('a, 'kind) linsolv_callbacks =
  | NoCallbacks

  | DenseCallback of dense_jac_fn
  | BandCallback  of band_jac_fn
  | SpilsCallback of 'a spils_callbacks

  | AlternateCallback of ('a, 'kind) alternate_linsolv

  | BBDCallback of 'a Bbd.callbacks

and ('a, 'k) session = {
  kinsol    : kin_mem;
  backref   : c_weak_ref;
  err_file  : kin_file;
  info_file : kin_file;

  mutable neqs       : int;    (* only valid for 'kind = serial *)
  mutable exn_temp   : exn option;

  mutable sysfn      : 'a -> 'a -> unit;
  mutable errh       : Sundials.error_details -> unit;
  mutable infoh      : Sundials.error_details -> unit;

  mutable ls_callbacks : ('a, 'k) linsolv_callbacks;
}

and ('data, 'kind) alternate_linsolv =
  {
    linit  : (('data, 'kind) session -> bool) option;
    lsetup : (('data, 'kind) session -> unit) option;
    lsolve : ('data, 'kind) session -> 'data -> 'data -> float;
  }


type ('data, 'kind) linear_solver = ('data, 'kind) session
                                        -> ('data, 'kind) nvector option -> unit

let read_weak_ref x : ('a, 'k) session =
  match Weak.get x 0 with
  | Some y -> y
  | None -> raise (Failure "Internal error: weak reference is dead")

let adjust_retcode = fun session check_recoverable f x ->
  try f x; 0
  with
  | Sundials.RecoverableFailure _ when check_recoverable -> 1
  | e -> (session.exn_temp <- Some e; -1)

let adjust_retcode_and_float = fun session f x ->
  try (f x, 0)
  with
  | Sundials.RecoverableFailure _ -> (0.0, 1)
  | e -> (session.exn_temp <- Some e; (0.0, -1))

