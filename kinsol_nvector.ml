(***********************************************************************)
(*                                                                     *)
(*               OCaml interface to (serial) Sundials                  *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

include Kinsol (* Need to register shared exceptions *)

type 'a nvector = 'a Nvector.nvector

type 'a single_tmp = 'a
type 'a double_tmp = 'a * 'a

type ('t, 'a) jacobian_arg =
  {
    jac_u   : 'a nvector;
    jac_fu  : 'a nvector;
    jac_tmp : 't
  }

type 'a linear_solver =
  | Spgmr of int option * int option * 'a spils_callbacks
  | Spbcg of int option * 'a spils_callbacks
  | Sptfqmr of int option * 'a spils_callbacks

and 'a spils_callbacks =
  {
    prec_solve_fn : (('a single_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                      -> 'a nvector -> unit) option;

    prec_setup_fn : (('a double_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                     -> unit) option;

    jac_times_vec_fn : ('a -> 'a -> 'a -> bool -> bool) option;
  }

and 'a prec_solve_arg = { uscale : 'a nvector; fscale : 'a nvector; }

let spils_no_precond = {
  prec_solve_fn = None;
  prec_setup_fn = None;
  jac_times_vec_fn = None;
}

type print_level =
  | NoInformation     (* 0 *)
  | ShowScaledNorms   (* 1 *)
  | ShowScaledDFNorm  (* 2 *)
  | ShowGlobalValues  (* 3 *)

type eta_params = {
  egamma : float option;
  ealpha : float option;
}

type eta_choice =
  | EtaChoice1                  (* KIN_ETACHOICE1 *)
  | EtaChoice2 of eta_params    (* KIN_ETACHOICE2 *)
  | EtaConstant of float option (* KIN_ETACONSTANT *)

type kin_mem
type kin_file
type c_weak_ref

type 'a session = {
  kinsol    : kin_mem;
  backref   : c_weak_ref;
  err_file  : kin_file;
  info_file : kin_file;

  mutable exn_temp : exn option;

  mutable sysfn      : 'a -> 'a -> unit;
  mutable errh       : Sundials.error_details -> unit;
  mutable infoh      : Sundials.error_details -> unit;
  mutable presetupfn : ('a double_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                       -> unit;
  mutable presolvefn : ('a single_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                       -> 'a nvector -> unit;
  mutable jactimesfn : 'a -> 'a -> 'a -> bool -> bool;
}

(* interface *)

let read_weak_ref x : 'a session =
  match Weak.get x 0 with
  | Some y -> y
  | None -> raise (Failure "Internal error: weak reference is dead")

let adjust_retcode = fun session check_recoverable f x ->
  try f x; 0
  with
  | Sundials.RecoverableFailure when check_recoverable -> 1
  | e -> (session.exn_temp <- Some e; -1)

let call_sysfn session u fval =
  let session = read_weak_ref session in
  adjust_retcode session true (session.sysfn u) fval

let call_errh session details =
  let session = read_weak_ref session in
  try session.errh details
  with e ->
    prerr_endline ("Warning: error handler function raised an exception.  " ^
                   "This exception will not be propagated.")

let call_infoh session details =
  let session = read_weak_ref session in
  try session.infoh details
  with e ->
    prerr_endline ("Warning: error handler function raised an exception.  " ^
                   "This exception will not be propagated.")

let call_presolvefn session jac ps =
  let session = read_weak_ref session in
  adjust_retcode session true (session.presolvefn jac) ps

let call_presetupfn session jac ps =
  let session = read_weak_ref session in
  adjust_retcode session true (session.presetupfn jac) ps

let call_jactimesfn session v jv u new_uu =
  let session = read_weak_ref session in
  try (session.jactimesfn v jv u new_uu, 0)
  with
  | Sundials.RecoverableFailure -> (false, 1)
  | e -> (session.exn_temp <- Some e; (false, -1))

let _ =
  Callback.register "c_nvec_kinsol_call_sysfn"         call_sysfn;
  Callback.register "c_nvec_kinsol_call_errh"          call_errh;
  Callback.register "c_nvec_kinsol_call_infoh"         call_infoh;
  Callback.register "c_nvec_kinsol_call_presolvefn"    call_presolvefn;
  Callback.register "c_nvec_kinsol_call_presetupfn"    call_presetupfn;
  Callback.register "c_nvec_kinsol_call_jactimesfn"    call_jactimesfn

external session_finalize : 'a session -> unit
    = "c_kinsol_session_finalize"

let shouldn't_be_called fcn =
  failwith ("internal error in sundials: " ^ fcn ^ " is called")
let dummy_dense_jac _ _ = shouldn't_be_called "dummy_dense_jac"
let dummy_band_jac _ _ _ _ = shouldn't_be_called "dummy_band_jac"
let dummy_prec_setup _ _ = shouldn't_be_called "dummy_prec_setup"
let dummy_prec_solve _ _ _ = shouldn't_be_called "dummy_prec_solve"
let dummy_jac_times_vec _ _ _ = shouldn't_be_called "dummy_jac_times_vec"

module Spils =
  struct
    type 'a solve_arg = 'a prec_solve_arg =
      {
        uscale : 'a nvector;
        fscale : 'a nvector;
      }

    (* TODO: CVODE/IDA versions to allow optional setup function. *)
    (* TODO: Eliminate this function and just provide reinit? *)
    external c_set_preconditioner : 'a session -> bool -> unit
        = "c_nvec_kinsol_spils_set_preconditioner"

    let set_preconditioner s fpresetupfn fpresolvefn =
      let withprec =
        match fpresetupfn with
        | None -> (s.presetupfn <- dummy_prec_setup; false)
        | Some fn -> (s.presetupfn <- fn; true)
      in
      s.presolvefn <- fpresolvefn;
      c_set_preconditioner s withprec

    (* TODO: Eliminate this function and just provide reinit? *)
    external c_set_jac_times_vec_fn : 'a session -> unit
        = "c_nvec_kinsol_spils_set_jac_times_vec_fn"

    let set_jac_times_vec_fn s fjactimesfn =
      s.jactimesfn <- fjactimesfn;
      c_set_jac_times_vec_fn s

    external c_clear_jac_times_vec_fn : 'a session -> unit
        = "c_nvec_cvode_clear_jac_times_vec_fn"

    let clear_jac_times_vec_fn s =
      s.jactimesfn <- dummy_jac_times_vec;
      c_clear_jac_times_vec_fn s

    external c_set_max_restarts     : 'a session -> int -> unit
        = "c_kinsol_spils_set_max_restarts"

    external get_work_space       : 'a session -> int * int
        = "c_kinsol_spils_get_work_space"

    external get_num_lin_iters    : 'a session -> int
        = "c_kinsol_spils_get_num_lin_iters"

    external get_num_conv_fails   : 'a session -> int
        = "c_kinsol_spils_get_num_conv_fails"

    external get_num_prec_evals   : 'a session -> int
        = "c_kinsol_spils_get_num_prec_evals"

    external get_num_prec_solves  : 'a session -> int
        = "c_kinsol_spils_get_num_prec_solves"

    external get_num_jtimes_evals : 'a session -> int
        = "c_kinsol_spils_get_num_jtimes_evals"

    external get_num_func_evals    : 'a session -> int
        = "c_kinsol_spils_get_num_func_evals"
  end

external set_error_file : 'a session -> string -> bool -> unit
    = "c_kinsol_set_error_file"

external c_set_err_handler_fn : 'a session -> unit
    = "c_nvec_kinsol_set_err_handler_fn"

let set_err_handler_fn s ferrh =
  s.errh <- ferrh;
  c_set_err_handler_fn s

external c_clear_err_handler_fn : 'a session -> unit
    = "c_nvec_kinsol_clear_err_handler_fn"

let clear_err_handler_fn s =
  s.errh <- (fun _ -> ());
  c_clear_err_handler_fn s

external set_info_file : 'a session -> string -> bool -> unit
    = "c_kinsol_set_info_file"

external c_set_info_handler_fn : 'a session -> unit
    = "c_nvec_kinsol_set_info_handler_fn"

let set_info_handler_fn s finfoh =
  s.infoh <- finfoh;
  c_set_info_handler_fn s

external c_clear_info_handler_fn : 'a session -> unit
    = "c_nvec_kinsol_clear_info_handler_fn"

let clear_info_handler_fn s =
  s.infoh <- (fun _ -> ());
  c_clear_info_handler_fn s

external set_print_level : 'a session -> print_level -> unit
    = "c_kinsol_set_print_level"

external set_num_max_iters : 'a session -> int -> unit
    = "c_kinsol_set_num_max_iters"

external set_no_init_setup : 'a session -> bool -> unit
    = "c_kinsol_set_no_init_setup"

let int_default = function None -> 0 | Some v -> v

external c_set_max_setup_calls : 'a session -> int -> unit
    = "c_kinsol_set_max_setup_calls"

let set_max_setup_calls s msbset =
  c_set_max_setup_calls s (match msbset with None -> 0 | Some i -> i)

external c_set_eta_form : 'a session -> eta_choice -> unit
    = "c_kinsol_set_eta_form"

external c_set_eta_const_value : 'a session -> float -> unit
    = "c_kinsol_set_eta_form"

external c_set_eta_params : 'a session -> float -> float -> unit
    = "c_kinsol_set_eta_form"

let float_default = function None -> 0.0 | Some v -> v

let set_eta_choice s etachoice =
  match etachoice with
  | EtaChoice1 -> ()
  | EtaChoice2 { egamma; ealpha } ->
      c_set_eta_params s (float_default egamma) (float_default ealpha)
  | EtaConstant eta ->
      c_set_eta_const_value s (float_default eta);
  c_set_eta_form s etachoice

external c_set_res_mon_const_value : 'a session -> float -> unit
    = "c_kinsol_set_res_mon_const_value"

let set_res_mon_const_value s omegaconst =
  c_set_res_mon_const_value s (float_default omegaconst)

external c_set_res_mon_params : 'a session -> float -> float -> unit
    = "c_kinsol_set_res_mon_params"

let set_res_mon_params s omegamin omegamax =
  c_set_res_mon_params s (float_default omegamin) (float_default omegamax)

external set_no_min_eps : 'a session -> bool -> unit
    = "c_kinsol_set_no_min_eps"

external c_set_max_newton_step : 'a session -> float -> unit
    = "c_kinsol_set_max_newton_step"

let set_max_newton_step s mxnewtstep =
  c_set_max_newton_step s (float_default mxnewtstep)

external c_set_max_beta_fails : 'a session -> float -> unit
    = "c_kinsol_set_max_beta_fails"

let set_max_beta_fails s mxnbcf =
  c_set_max_beta_fails s (float_default mxnbcf)

external c_set_rel_err_func : 'a session -> float -> unit
    = "c_kinsol_set_rel_err_func"

let set_rel_err_func s relfunc =
  c_set_rel_err_func s (float_default relfunc)

external c_set_func_norm_tol : 'a session -> float -> unit
    = "c_kinsol_set_func_norm_tol"

let set_func_norm_tol s fnormtol =
  c_set_func_norm_tol s (float_default fnormtol)

external c_set_scaled_step_tol : 'a session -> float -> unit
    = "c_kinsol_set_scaled_step_tol"

let set_scaled_step_tol s scsteptol =
  c_set_scaled_step_tol s (float_default scsteptol)

external set_constraints : 'a session -> 'a -> unit
    = "c_nvec_kinsol_set_constraints"

let set_sys_func s fsys =
  s.sysfn <- fsys

external get_work_space : 'a session -> int * int
    = "c_kinsol_get_work_space"

external get_num_func_evals : 'a session -> int
    = "c_kinsol_get_num_func_evals"

external get_num_nonlin_solv_iters : 'a session -> int
    = "c_kinsol_get_num_nonlin_solv_iters"

external get_num_beta_cond_fails : 'a session -> int
    = "c_kinsol_get_num_beta_cond_fails"

external get_num_backtrack_ops : 'a session -> int
    = "c_kinsol_get_num_backtrack_ops"

external get_func_norm : 'a session -> float
    = "c_kinsol_get_func_norm"

external get_step_length : 'a session -> float
    = "c_kinsol_get_step_length"

external c_spils_spgmr : 'a session -> int -> unit
  = "c_kinsol_spils_spgmr"

external c_spils_spbcg : 'a session -> int -> unit
  = "c_kinsol_spils_spbcg"

external c_spils_sptfqmr : 'a session -> int -> unit
  = "c_kinsol_spils_sptfqmr"

let set_spils_callbacks s {prec_solve_fn; prec_setup_fn; jac_times_vec_fn} =
  match (prec_solve_fn, prec_setup_fn) with
  | None, _ -> ()
  | Some solve, osetup -> Spils.set_preconditioner s osetup solve;
  match jac_times_vec_fn with
  | None -> ()
  | Some jtimes -> Spils.set_jac_times_vec_fn s jtimes

let set_linear_solver s lsolver =
  match lsolver with
  | Spgmr (maxl, omaxrs, callbacks) ->
      c_spils_spgmr s (int_default maxl);
      (match omaxrs with
       | None -> ()
       | Some maxrs -> Spils.c_set_max_restarts s maxrs);
      set_spils_callbacks s callbacks

  | Spbcg (maxl, callbacks) ->
      c_spils_spbcg s (int_default maxl);
      set_spils_callbacks s callbacks

  | Sptfqmr (maxl, callbacks) ->
      c_spils_sptfqmr s (int_default maxl);
      set_spils_callbacks s callbacks

external c_init
    : 'a session Weak.t -> 'a nvector
      -> (kin_mem * c_weak_ref * kin_file * kin_file)
    = "c_nvec_kinsol_init"

let init lsolver f u0 =
  let weakref = Weak.create 1 in
  let kin_mem, backref, err_file, info_file = c_init weakref u0
  in
  let session = {
          kinsol     = kin_mem;
          backref    = backref;
          err_file   = err_file;
          info_file  = info_file;

          exn_temp   = None;

          sysfn      = f;
          errh       = (fun _ -> ());
          infoh      = (fun _ -> ());
          presetupfn = dummy_prec_setup;
          presolvefn = dummy_prec_solve;
          jactimesfn = dummy_jac_times_vec;
        } in
  Gc.finalise session_finalize session;
  Weak.set weakref 0 (Some session);
  set_linear_solver session lsolver;
  session

type result =
  | Success           (** KIN_SUCCESS *)
  | InitialGuessOK    (** KIN_INITIAL_GUESS_OK *)
  | StoppedOnStepTol  (** KIN_STEP_LT_STPTOL *)

external solve : 'a session -> 'a nvector -> bool -> 'a nvector
                  -> 'a nvector -> result
    = "c_nvec_kinsol_solve"

