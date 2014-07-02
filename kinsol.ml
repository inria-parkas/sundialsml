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

let _ =
  List.iter (fun (nm, ex) -> Callback.register_exception nm ex)
  [
    ("kinsol_RecoverableFailure",             Sundials.RecoverableFailure);

    ("kinsol_IllInput",                       IllInput);
    ("kinsol_LineSearchNonConvergence",       LineSearchNonConvergence);
    ("kinsol_MaxIterationsReached",           MaxIterationsReached);
    ("kinsol_MaxNewtonStepExceeded",          MaxNewtonStepExceeded);
    ("kinsol_LineSearchBetaConditionFailure", LineSearchBetaConditionFailure);
    ("kinsol_LinearSolverNoRecovery",         LinearSolverNoRecovery);
    ("kinsol_LinearSolverInitFailure",        LinearSolverInitFailure);
    ("kinsol_LinearSetupFailure",             LinearSetupFailure);
    ("kinsol_LinearSolverFailure",            LinearSolverFailure);
    ("kinsol_SystemFunctionFailure",          SystemFunctionFailure);
    ("kinsol_FirstSystemFunctionFailure",     FirstSystemFunctionFailure);
    ("kinsol_RepeatedSystemFunctionFailure",  RepeatedSystemFunctionFailure);
  ]

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

type 'a prec_solve_arg = { uscale : 'a; fscale : 'a; }

type 'a spils_callbacks =
  {
    prec_solve_fn : (('a single_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                      -> 'a -> unit) option;

    prec_setup_fn : (('a double_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                     -> unit) option;

    jac_times_vec_fn : ('a -> 'a -> 'a -> bool -> bool) option;
  }

type ('a, 'kind) linsolv_callbacks =
  | NoCallbacks

  | DenseCallback of dense_jac_fn
  | BandCallback  of band_jac_fn
  | SpilsCallback of 'a spils_callbacks

type ('a, 'k) session = {
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

type serial_session = (real_array, Nvector_serial.kind) session

type ('data, 'kind) linear_solver =
  ('data, 'kind) session -> (('data, 'kind) nvector) option -> unit
type serial_linear_solver = (real_array, Nvector_serial.kind) linear_solver

(* interface *)

let read_weak_ref x : ('a, 'k) session =
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

let call_jacfn session jac j =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | DenseCallback f -> adjust_retcode session true (f jac) j
  | _ -> assert false

let call_bandjacfn session range jac j =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | BandCallback f -> adjust_retcode session true (f range jac) j
  | _ -> assert false

let call_precsolvefn session jac ps u =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | SpilsCallback { prec_solve_fn = Some f } ->
      adjust_retcode session true (f jac ps) u
  | _ -> assert false

let call_precsetupfn session jac ps =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | SpilsCallback { prec_setup_fn = Some f } ->
      adjust_retcode session true (f jac) ps
  | _ -> assert false

let call_jactimesfn session v jv u new_uu =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | SpilsCallback { jac_times_vec_fn = Some f } ->
      (try (f v jv u new_uu, 0) with
       | Sundials.RecoverableFailure -> (false, 1)
       | e -> (session.exn_temp <- Some e; (false, -1)))
  | _ -> assert false

let _ =
  Callback.register "c_kinsol_call_sysfn"         call_sysfn;
  Callback.register "c_kinsol_call_errh"          call_errh;
  Callback.register "c_kinsol_call_infoh"         call_infoh;
  Callback.register "c_kinsol_call_jacfn"         call_jacfn;
  Callback.register "c_kinsol_call_bandjacfn"     call_bandjacfn;
  Callback.register "c_kinsol_call_precsolvefn"   call_precsolvefn;
  Callback.register "c_kinsol_call_precsetupfn"   call_precsetupfn;
  Callback.register "c_kinsol_call_jactimesfn"    call_jactimesfn

external session_finalize : ('a, 'k) session -> unit
    = "c_kinsol_session_finalize"

let int_default = function None -> 0 | Some v -> v

module Dls =
  struct
    type dense_jac_fn = (real_array double_tmp, real_array) jacobian_arg
                                                -> Dls.DenseMatrix.t -> unit

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
       | Some nv -> s.neqs <- Sundials.RealArray.length (Sundials.unvec nv)
       | None -> ());
      c_dls_dense s (fo <> None);
      s.ls_callbacks <- match fo with
                        | None -> NoCallbacks
                        | Some f -> DenseCallback f

    let lapack_dense fo s onv =
      (match onv with
       | Some nv -> s.neqs <- Sundials.RealArray.length (Sundials.unvec nv)
       | None -> ());
      c_dls_lapack_dense s (fo <> None);
      s.ls_callbacks <- match fo with
                        | None -> NoCallbacks
                        | Some f -> DenseCallback f

    type band_jac_fn = bandrange
                        -> (real_array double_tmp, real_array) jacobian_arg
                        -> Dls.BandMatrix.t -> unit

    let band { mupper; mlower } fo s onv =
      (match onv with
       | Some nv -> s.neqs <- Sundials.RealArray.length (Sundials.unvec nv)
       | None -> ());
      c_dls_band s mupper mlower (fo <> None);
      s.ls_callbacks <- match fo with
                        | None -> NoCallbacks
                        | Some f -> BandCallback f

    let lapack_band { mupper; mlower } fo s onv =
      (match onv with
       | Some nv -> s.neqs <- Sundials.RealArray.length (Sundials.unvec nv)
       | None -> ());
      c_dls_lapack_band s mupper mlower (fo <> None);
      s.ls_callbacks <- match fo with
                        | None -> NoCallbacks
                        | Some f -> BandCallback f

    external set_dense_jac_fn : serial_session -> unit
        = "c_kinsol_dls_set_dense_jac_fn"

    let set_dense_jac_fn s fjacfn =
      s.ls_callbacks <- DenseCallback fjacfn;
      set_dense_jac_fn s

    external clear_dense_jac_fn : serial_session -> unit
        = "c_kinsol_dls_clear_dense_jac_fn"

    let clear_dense_jac_fn s =
      match s.ls_callbacks with
      | DenseCallback _ -> (s.ls_callbacks <- NoCallbacks;
                            clear_dense_jac_fn s)
      | _ -> failwith "dense linear solver not in use"

    external set_band_jac_fn : serial_session -> unit
        = "c_kinsol_dls_set_band_jac_fn"

    let set_band_jac_fn s fbandjacfn =
      s.ls_callbacks <- BandCallback fbandjacfn;
      set_band_jac_fn s

    external clear_band_jac_fn : serial_session -> unit
        = "c_kinsol_dls_clear_band_jac_fn"

    let clear_band_jac_fn s =
      match s.ls_callbacks with
      | BandCallback _ -> (s.ls_callbacks <- NoCallbacks;
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
    type 'a solve_arg = 'a prec_solve_arg =
      {
        uscale : 'a;
        fscale : 'a;
      }

    type 'a callbacks = 'a spils_callbacks =
      {
        prec_solve_fn : (('a single_tmp, 'a) jacobian_arg -> 'a solve_arg
                          -> 'a -> unit) option;

        prec_setup_fn : (('a double_tmp, 'a) jacobian_arg -> 'a solve_arg
                         -> unit) option;

        jac_times_vec_fn : ('a -> 'a -> 'a -> bool -> bool) option;
      }

    let spils_no_precond = {
      prec_solve_fn = None;
      prec_setup_fn = None;
      jac_times_vec_fn = None;
    }

    external c_set_max_restarts     : ('a, 'k) session -> int -> unit
        = "c_kinsol_spils_set_max_restarts"

    external c_set_preconditioner : ('a, 'k) session -> bool -> unit
        = "c_kinsol_spils_set_preconditioner"

    let set_preconditioner s fprecsetupfn fprecsolvefn =
      (match s.ls_callbacks with
       | SpilsCallback cbs ->
           s.ls_callbacks <- SpilsCallback { cbs with
                               prec_setup_fn = fprecsetupfn;
                               prec_solve_fn = Some fprecsolvefn }
       | _ -> failwith "spils solver not in use");
      c_set_preconditioner s (fprecsetupfn <> None)

    external c_set_jac_times_vec_fn : ('a, 'k) session -> unit
        = "c_kinsol_spils_set_jac_times_vec_fn"

    let set_jac_times_vec_fn s fjactimesfn =
      (match s.ls_callbacks with
       | SpilsCallback cbs ->
           s.ls_callbacks <- SpilsCallback { cbs with
                               jac_times_vec_fn = Some fjactimesfn }
       | _ -> failwith "spils solver not in use");
      c_set_jac_times_vec_fn s

    let set_callbacks s ({prec_solve_fn; prec_setup_fn; jac_times_vec_fn} as cb)
      = (match (prec_solve_fn, prec_setup_fn) with
         | None, _ -> ()
         | Some solve, osetup -> c_set_preconditioner s (osetup <> None));
        (match jac_times_vec_fn with
         | None -> ()
          | Some jtimes -> c_set_jac_times_vec_fn s);
        s.ls_callbacks <- SpilsCallback cb

    external c_spils_spgmr : ('a, 'k) session -> int -> unit
      = "c_kinsol_spils_spgmr"

    let spgmr maxl omaxrs callbacks s _ =
      c_spils_spgmr s (int_default maxl);
      (match omaxrs with
       | None -> ()
       | Some maxrs -> c_set_max_restarts s maxrs);
      set_callbacks s callbacks

    external c_spils_spbcg : ('a, 'k) session -> int -> unit
      = "c_kinsol_spils_spbcg"

    let spbcg maxl callbacks s _ =
      c_spils_spbcg s (int_default maxl);
      set_callbacks s callbacks

    external c_spils_sptfqmr : ('a, 'k) session -> int -> unit
      = "c_kinsol_spils_sptfqmr"

    let sptfqmr maxl callbacks s _ =
      c_spils_sptfqmr s (int_default maxl);
      set_callbacks s callbacks

    external c_clear_jac_times_vec_fn : ('a, 'k) session -> unit
        = "c_cvode_clear_jac_times_vec_fn"

    let clear_jac_times_vec_fn s =
      (match s.ls_callbacks with
       | SpilsCallback cbs ->
           s.ls_callbacks <- SpilsCallback { cbs with jac_times_vec_fn = None }
       | _ -> failwith "spils solver not in use");
      c_clear_jac_times_vec_fn s

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
  s.errh <- (fun _ -> ());
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
  s.infoh <- (fun _ -> ());
  c_clear_info_handler_fn s

external set_print_level : ('a, 'k) session -> print_level -> unit
    = "c_kinsol_set_print_level"

external set_num_max_iters : ('a, 'k) session -> int -> unit
    = "c_kinsol_set_num_max_iters"

external set_no_init_setup : ('a, 'k) session -> bool -> unit
    = "c_kinsol_set_no_init_setup"

external set_no_res_mon : ('a, 'k) session -> bool -> unit
    = "c_kinsol_set_no_res_mon"

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

let set_res_mon_params s omegamin omegamax =
  c_set_res_mon_params s (float_default omegamin) (float_default omegamax)

external set_no_min_eps : ('a, 'k) session -> bool -> unit
    = "c_kinsol_set_no_min_eps"

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
          errh         = (fun _ -> ());
          infoh        = (fun _ -> ());

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

