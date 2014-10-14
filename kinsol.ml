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

include Kinsol_impl
open Sundials

exception IllInput                       (* KIN_ILL_INPUT *)
exception LineSearchNonConvergence       (* KIN_LINESEARCH_NONCONV *)
exception MaxIterationsReached           (* KIN_MAXITER_REACHED *)
exception MaxNewtonStepExceeded          (* KIN_MXNEWT_5X_EXCEEDED *)
exception LineSearchBetaConditionFailure (* KIN_LINESEARCH_BCFAIL *)
exception LinearSolverNoRecovery         (* KIN_LINSOLV_NO_RECOVERY *)
exception LinearSolverInitFailure        (* KIN_LINIT_FAIL *)
exception LinearSetupFailure             (* KIN_LSETUP_FAIL *)
exception LinearSolverFailure            (* KIN_LSOLVE_FAIL *)
exception SystemFunctionFailure          (* KIN_SYSFUNC_FAIL *)
exception FirstSystemFunctionFailure     (* KIN_FIRST_SYSFUNC_FAIL *)
exception RepeatedSystemFunctionFailure  (* KIN_REPTD_SYSFUNC_ERR *)

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

type serial_linear_solver = (RealArray.t, Nvector_serial.kind) linear_solver
type 'a double = 'a * 'a

(* interface *)

let int_default = function None -> 0 | Some v -> v

module Dls =
  struct
    include DlsTypes

    external c_dls_dense : serial_session -> bool -> unit
      = "c_kinsol_dls_dense"

    external c_dls_lapack_dense : serial_session -> bool -> unit
      = "c_kinsol_dls_lapack_dense"

    external c_dls_band : serial_session -> int -> int -> bool -> unit
      = "c_kinsol_dls_band"

    external c_dls_lapack_band : serial_session -> int -> int -> bool -> unit
      = "c_kinsol_dls_lapack_band"

    let dense fo s onv =
      (match onv with
       | Some nv -> s.neqs <- Sundials.RealArray.length (Nvector.unwrap nv)
       | None -> ());
      c_dls_dense s (fo <> None);
      s.ls_callbacks <- match fo with
                        | None -> NoCallbacks
                        | Some f -> DenseCallback { jacfn = f; dmat = None }

    let lapack_dense fo s onv =
      (match onv with
       | Some nv -> s.neqs <- Sundials.RealArray.length (Nvector.unwrap nv)
       | None -> ());
      c_dls_lapack_dense s (fo <> None);
      s.ls_callbacks <- match fo with
                        | None -> NoCallbacks
                        | Some f -> DenseCallback { jacfn = f; dmat = None }

    let band { mupper; mlower } fo s onv =
      (match onv with
       | Some nv -> s.neqs <- Sundials.RealArray.length (Nvector.unwrap nv)
       | None -> ());
      c_dls_band s mupper mlower (fo <> None);
      s.ls_callbacks <- match fo with
                        | None -> NoCallbacks
                        | Some f -> BandCallback { bjacfn = f; bmat = None }

    let lapack_band { mupper; mlower } fo s onv =
      (match onv with
       | Some nv -> s.neqs <- Sundials.RealArray.length (Nvector.unwrap nv)
       | None -> ());
      c_dls_lapack_band s mupper mlower (fo <> None);
      s.ls_callbacks <- match fo with
                        | None -> NoCallbacks
                        | Some f -> BandCallback { bjacfn = f; bmat = None }

    let invalidate_callback (type d) (type k) (session : (d, k) session) =
      match session.ls_callbacks with
      | DenseCallback ({ dmat = Some d } as cb) ->
          Dls.DenseMatrix.invalidate d;
          cb.dmat <- None
      | BandCallback  ({ bmat = Some d } as cb) ->
          Dls.BandMatrix.invalidate d;
          cb.bmat <- None
      | _ -> ()

    external set_dense_jac_fn : serial_session -> unit
        = "c_kinsol_dls_set_dense_jac_fn"

    let set_dense_jac_fn s fjacfn =
      invalidate_callback s;
      s.ls_callbacks <- DenseCallback { jacfn = fjacfn; dmat = None };
      set_dense_jac_fn s

    external clear_dense_jac_fn : serial_session -> unit
        = "c_kinsol_dls_clear_dense_jac_fn"

    let clear_dense_jac_fn s =
      match s.ls_callbacks with
      | DenseCallback _ -> (invalidate_callback s;
                            s.ls_callbacks <- NoCallbacks;
                            clear_dense_jac_fn s)
      | _ -> failwith "dense linear solver not in use"

    external set_band_jac_fn : serial_session -> unit
        = "c_kinsol_dls_set_band_jac_fn"

    let set_band_jac_fn s fbandjacfn =
      invalidate_callback s;
      s.ls_callbacks <- BandCallback { bjacfn = fbandjacfn; bmat = None };
      set_band_jac_fn s

    external clear_band_jac_fn : serial_session -> unit
        = "c_kinsol_dls_clear_band_jac_fn"

    let clear_band_jac_fn s =
      match s.ls_callbacks with
      | BandCallback _ -> (invalidate_callback s;
                           s.ls_callbacks <- NoCallbacks;
                           clear_band_jac_fn s)
      | _ -> failwith "dense linear solver not in use"

    external get_work_space : serial_session -> int * int
        = "c_kinsol_dls_get_work_space"

    external get_num_jac_evals : serial_session -> int
        = "c_kinsol_dls_get_num_jac_evals"

    external get_num_func_evals : serial_session -> int
        = "c_kinsol_dls_get_num_func_evals"
  end

module Spils =
  struct
    include SpilsTypes

    external c_spgmr : ('a, 'k) session -> int -> unit
      = "c_kinsol_spils_spgmr"

    external c_spbcg : ('a, 'k) session -> int -> unit
      = "c_kinsol_spils_spbcg"

    external c_sptfqmr : ('a, 'k) session -> int -> unit
      = "c_kinsol_spils_sptfqmr"

    external c_set_max_restarts     : ('a, 'k) session -> int -> unit
        = "c_kinsol_spils_set_max_restarts"

    external c_set_preconditioner : ('a, 'k) session -> bool -> bool -> unit
        = "c_kinsol_spils_set_preconditioner"

    external c_set_jac_times_vec_fn : ('a, 'k) session -> bool -> unit
        = "c_kinsol_spils_set_jac_times_vec_fn"

    let init_preconditioner solve setup jac_times session nv =
      c_set_preconditioner session (solve <> None) (setup <> None);
      c_set_jac_times_vec_fn session (jac_times <> None);
      session.ls_callbacks <-
        SpilsCallback { prec_solve_fn = solve;
                        prec_setup_fn = setup;
                        jac_times_vec_fn = jac_times }

    let prec_none = InternalPrecNone
    let prec_right ?setup ?solve ?jac_times_vec () =
      InternalPrecRight (init_preconditioner solve setup jac_times_vec)

    let init_spils init maxl prec session nv =
      match prec with
      | InternalPrecNone -> init session maxl
      | InternalPrecRight set_prec ->
        init session maxl;
        set_prec session nv

    let spgmr ?(maxl=0) ?(max_restarts=5) prec session nv =
      init_spils c_spgmr maxl prec session nv;
      (* Note: we can skip set_max_restarts only when initializing a
         fresh solver.  *)
      if max_restarts <> 5 then
        c_set_max_restarts session max_restarts

    let spbcg ?(maxl=0) prec session nv =
      init_spils c_spbcg maxl prec session nv

    let sptfqmr ?(maxl=0) prec session nv =
      init_spils c_sptfqmr maxl prec session nv

    let set_preconditioner s ?setup ?solve () =
      match s.ls_callbacks with
      | SpilsCallback cbs ->
        s.ls_callbacks <- SpilsCallback { cbs with
                                          prec_setup_fn = setup;
                                          prec_solve_fn = solve };
        c_set_preconditioner s (solve <> None) (setup <> None)
      | _ -> failwith "spils solver not in use"

    let set_jac_times_vec_fn s f =
      match s.ls_callbacks with
      | SpilsCallback cbs ->
        s.ls_callbacks <- SpilsCallback { cbs with jac_times_vec_fn = Some f };
        c_set_jac_times_vec_fn s true
      | _ -> failwith "spils solver not in use"

    let clear_jac_times_vec_fn s =
      match s.ls_callbacks with
      | SpilsCallback cbs ->
        s.ls_callbacks <- SpilsCallback { cbs with jac_times_vec_fn = None };
        c_set_jac_times_vec_fn s false
      | _ -> failwith "spils solver not in use"

    external get_work_space       : ('a, 'k) session -> int * int
        = "c_kinsol_spils_get_work_space"

    external get_num_lin_iters    : ('a, 'k) session -> int
        = "c_kinsol_spils_get_num_lin_iters"

    external get_num_conv_fails   : ('a, 'k) session -> int
        = "c_kinsol_spils_get_num_conv_fails"

    external get_num_prec_evals   : ('a, 'k) session -> int
        = "c_kinsol_spils_get_num_prec_evals"

    external get_num_prec_solves  : ('a, 'k) session -> int
        = "c_kinsol_spils_get_num_prec_solves"

    external get_num_jtimes_evals : ('a, 'k) session -> int
        = "c_kinsol_spils_get_num_jtimes_evals"

    external get_num_func_evals    : ('a, 'k) session -> int
        = "c_kinsol_spils_get_num_func_evals"
  end

module Alternate =
  struct
    include AlternateTypes

    external c_set_alternate
      : ('data, 'kind) session -> bool -> bool -> unit
      = "c_kinsol_set_alternate"

    external get_u_uscale : ('data, 'kind) session -> 'data * 'data
      = "c_kinsol_get_u_uscale"

    external get_f_fscale  : ('data, 'kind) session -> 'data * 'data
      = "c_kinsol_get_f_fscale"

    external set_sjpnorm   : ('data, 'kind) session -> float -> unit
      = "c_kinsol_set_sjpnorm"

    external set_sfdotjp   : ('data, 'kind) session -> float -> unit
      = "c_kinsol_set_sfdotjp"

    let make f s nv =
      let { linit; lsetup; lsolve } as cb = f s nv in
      c_set_alternate s (linit <> None) (lsetup <> None);
      s.ls_callbacks <- AlternateCallback cb

  end

external set_error_file : ('a, 'k) session -> string -> bool -> unit
    = "c_kinsol_set_error_file"

external c_set_err_handler_fn : ('a, 'k) session -> unit
    = "c_kinsol_set_err_handler_fn"

let set_err_handler_fn s ferrh =
  s.errh <- ferrh;
  c_set_err_handler_fn s

external c_clear_err_handler_fn : ('a, 'k) session -> unit
    = "c_kinsol_clear_err_handler_fn"

let clear_err_handler_fn s =
  s.errh <- dummy_errh;
  c_clear_err_handler_fn s

external set_info_file : ('a, 'k) session -> string -> bool -> unit
    = "c_kinsol_set_info_file"

external c_set_info_handler_fn : ('a, 'k) session -> unit
    = "c_kinsol_set_info_handler_fn"

let set_info_handler_fn s finfoh =
  s.infoh <- finfoh;
  c_set_info_handler_fn s

external c_clear_info_handler_fn : ('a, 'k) session -> unit
    = "c_kinsol_clear_info_handler_fn"

let clear_info_handler_fn s =
  s.infoh <- dummy_infoh;
  c_clear_info_handler_fn s

external set_print_level : ('a, 'k) session -> print_level -> unit
    = "c_kinsol_set_print_level"

external set_num_max_iters : ('a, 'k) session -> int -> unit
    = "c_kinsol_set_num_max_iters"

external c_set_no_init_setup : ('a, 'k) session -> bool -> unit
    = "c_kinsol_set_no_init_setup"

let set_no_init_setup s = c_set_no_init_setup s true
let set_init_setup s = c_set_no_init_setup s false

external c_set_no_res_mon : ('a, 'k) session -> bool -> unit
    = "c_kinsol_set_no_res_mon"

let set_no_res_mon s = c_set_no_res_mon s true
let set_res_mon s = c_set_no_res_mon s false

external c_set_max_sub_setup_calls : ('a, 'k) session -> int -> unit
    = "c_kinsol_set_max_sub_setup_calls"

let set_max_sub_setup_calls s msbset =
  c_set_max_sub_setup_calls s (int_default msbset)

external c_set_max_setup_calls : ('a, 'k) session -> int -> unit
    = "c_kinsol_set_max_setup_calls"

let set_max_setup_calls s msbset =
  c_set_max_setup_calls s (match msbset with None -> 0 | Some i -> i)

external c_set_eta_form : ('a, 'k) session -> eta_choice -> unit
    = "c_kinsol_set_eta_form"

external c_set_eta_const_value : ('a, 'k) session -> float -> unit
    = "c_kinsol_set_eta_form"

external c_set_eta_params : ('a, 'k) session -> float -> float -> unit
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

external c_set_res_mon_const_value : ('a, 'k) session -> float -> unit
    = "c_kinsol_set_res_mon_const_value"

let set_res_mon_const_value s omegaconst =
  c_set_res_mon_const_value s (float_default omegaconst)

external c_set_res_mon_params : ('a, 'k) session -> float -> float -> unit
    = "c_kinsol_set_res_mon_params"

let set_res_mon_params s ?omegamin ?omegamax () =
  c_set_res_mon_params s (float_default omegamin) (float_default omegamax)

external c_set_no_min_eps : ('a, 'k) session -> bool -> unit
    = "c_kinsol_set_no_min_eps"

let set_no_min_eps s = c_set_no_min_eps s true
let set_min_eps s = c_set_no_min_eps s false

external c_set_max_newton_step : ('a, 'k) session -> float -> unit
    = "c_kinsol_set_max_newton_step"

let set_max_newton_step s mxnewtstep =
  c_set_max_newton_step s (float_default mxnewtstep)

external c_set_max_beta_fails : ('a, 'k) session -> float -> unit
    = "c_kinsol_set_max_beta_fails"

let set_max_beta_fails s mxnbcf =
  c_set_max_beta_fails s (float_default mxnbcf)

external c_set_rel_err_func : ('a, 'k) session -> float -> unit
    = "c_kinsol_set_rel_err_func"

let set_rel_err_func s relfunc =
  c_set_rel_err_func s (float_default relfunc)

external c_set_func_norm_tol : ('a, 'k) session -> float -> unit
    = "c_kinsol_set_func_norm_tol"

let set_func_norm_tol s fnormtol =
  c_set_func_norm_tol s (float_default fnormtol)

external c_set_scaled_step_tol : ('a, 'k) session -> float -> unit
    = "c_kinsol_set_scaled_step_tol"

let set_scaled_step_tol s scsteptol =
  c_set_scaled_step_tol s (float_default scsteptol)

module Constraint =
  struct
    let unconstrained = 0.0
    let non_negative = 1.0
    let non_positive = -1.0
    let positive = 2.0
    let negative = -2.0

    type t =
    | Unconstrained
    | NonNegative
    | NonPositive
    | Positive
    | Negative

    let of_float = function
      | 0.0  -> Unconstrained
      | 1.0  -> NonNegative
      | -1.0 -> NonPositive
      | 2.0  -> Positive
      | -2.0 -> Negative
      | f -> raise (Invalid_argument
                      ("invalid constraint: " ^ string_of_float f))
    let to_float = function
      | Unconstrained -> 0.0
      | NonNegative   -> 1.0
      | NonPositive   -> -1.0
      | Positive      -> 2.0
      | Negative      -> -2.0

    let name_of_constraint = function
      | Unconstrained -> "Unconstrained"
      | NonNegative -> "NonNegative"
      | NonPositive -> "NonPositive"
      | Positive -> "Positive"
      | Negative -> "Negative"

    let name_of_float x = name_of_constraint (of_float x)

    let string_of_constraint = function
      | Unconstrained -> invalid_arg "unconstrained"
      | NonNegative -> ">= 0"
      | NonPositive -> "<= 0"
      | Positive -> "> 0"
      | Negative -> "< 0"

    let string_of_float x = string_of_constraint (of_float x)
  end

external set_constraints : ('a, 'k) session -> ('a, 'k) nvector -> unit
    = "c_kinsol_set_constraints"

let set_linear_solver s lin_solv = lin_solv s None

let set_sys_func s fsys =
  s.sysfn <- fsys

external get_work_space : ('a, 'k) session -> int * int
    = "c_kinsol_get_work_space"

external get_num_func_evals : ('a, 'k) session -> int
    = "c_kinsol_get_num_func_evals"

external get_num_nonlin_solv_iters : ('a, 'k) session -> int
    = "c_kinsol_get_num_nonlin_solv_iters"

external get_num_beta_cond_fails : ('a, 'k) session -> int
    = "c_kinsol_get_num_beta_cond_fails"

external get_num_backtrack_ops : ('a, 'k) session -> int
    = "c_kinsol_get_num_backtrack_ops"

external get_func_norm : ('a, 'k) session -> float
    = "c_kinsol_get_func_norm"

external get_step_length : ('a, 'k) session -> float
    = "c_kinsol_get_step_length"

external c_init
    : ('a, 'k) session Weak.t -> ('a, 'k) nvector
      -> (kin_mem * c_weak_ref * kin_file * kin_file)
    = "c_kinsol_init"

external c_session_finalize : ('a, 'k) session -> unit
    = "c_kinsol_session_finalize"

let session_finalize s =
  Dls.invalidate_callback s;
  c_session_finalize s

let init lsolver f u0 =
  let weakref = Weak.create 1 in
  let kin_mem, backref, err_file, info_file = c_init weakref u0
  in
  let session = {
          kinsol       = kin_mem;
          backref      = backref;
          err_file     = err_file;
          info_file    = info_file;

          exn_temp     = None;

          neqs         = 0;

          sysfn        = f;
          errh         = dummy_errh;
          infoh        = dummy_infoh;

          ls_callbacks = NoCallbacks
        } in
  Gc.finalise session_finalize session;
  Weak.set weakref 0 (Some session);
  lsolver session (Some u0);
  session

type result =
  | Success           (** KIN_SUCCESS *)
  | InitialGuessOK    (** KIN_INITIAL_GUESS_OK *)
  | StoppedOnStepTol  (** KIN_STEP_LT_STPTOL *)

external solve : ('a, 'k) session -> ('a, 'k) nvector -> bool -> ('a, 'k) nvector
                  -> ('a, 'k) nvector -> result
    = "c_kinsol_solve"


(* Callbacks *)

let call_errh session details =
  let session = read_weak_ref session in
  try session.errh details
  with e ->
    prerr_endline ("Warning: error handler function raised an exception.  " ^
                   "This exception will not be propagated: " ^
                   Printexc.to_string e)

let call_infoh session details =
  let session = read_weak_ref session in
  try session.infoh details
  with e ->
    prerr_endline ("Warning: info function raised an exception.  " ^
                   "This exception will not be propagated: " ^
                   Printexc.to_string e)

(* Let C code know about some of the values in this module.  *)
external c_init_module : 'fcns -> exn array -> unit =
  "c_kinsol_init_module"

let _ =
  c_init_module
    (* Functions must be listed in the same order as
       callback_index in kinsol_ml.c.  *)
    (call_errh, call_infoh)

    (* Exceptions must be listed in the same order as
       kinsol_exn_index.  *)
    [|
      IllInput;
      LineSearchNonConvergence;
      MaxIterationsReached;
      MaxNewtonStepExceeded;
      LineSearchBetaConditionFailure;
      LinearSolverNoRecovery;
      LinearSolverInitFailure;
      LinearSetupFailure;
      LinearSolverFailure;
      SystemFunctionFailure;
      FirstSystemFunctionFailure;
      RepeatedSystemFunctionFailure;
    |]
