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
include Kinsol_impl

(* "Simulate" Linear Solvers in Sundials < 3.0.0 *)
let in_compat_mode =
  match Config.sundials_version with
  | 2,_,_ -> true
  | _ -> false

exception IllInput                       (* KIN_ILL_INPUT *)
exception LineSearchNonConvergence       (* KIN_LINESEARCH_NONCONV *)
exception MaxIterationsReached           (* KIN_MAXITER_REACHED *)
exception MaxNewtonStepExceeded          (* KIN_MXNEWT_5X_EXCEEDED *)
exception LineSearchBetaConditionFailure (* KIN_LINESEARCH_BCFAIL *)
exception LinearSolverNoRecovery         (* KIN_LINSOLV_NO_RECOVERY *)
exception LinearSolverInitFailure        (* KIN_LINIT_FAIL *)
exception LinearSetupFailure of exn option (* KIN_LSETUP_FAIL *)
exception LinearSolveFailure of exn option (* KIN_LSOLVE_FAIL *)
exception SystemFunctionFailure          (* KIN_SYSFUNC_FAIL *)
exception FirstSystemFunctionFailure     (* KIN_FIRST_SYSFUNC_FAIL *)
exception RepeatedSystemFunctionFailure  (* KIN_REPTD_SYSFUNC_ERR *)
exception VectorOpErr                    (* KIN_VECTOROP_ERR *)
exception MissingLinearSolver

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

type 'k serial_session_linear_solver =
  (RealArray.t, 'k) session_linear_solver
  constraint 'k = [>Nvector_serial.kind]

(* interface *)

let int_default = function None -> 0 | Some v -> v

module Dls = struct (* {{{ *)
  include DirectTypes
  include LinearSolver.Direct

  (* Sundials < 3.0.0 *)
  external c_dls_dense : 'k serial_session -> bool -> unit
    = "sunml_kinsol_dls_dense"

  (* Sundials < 3.0.0 *)
  external c_dls_lapack_dense : 'k serial_session -> bool -> unit
    = "sunml_kinsol_dls_lapack_dense"

  (* Sundials < 3.0.0 *)
  external c_dls_band : 'k serial_session -> int -> int -> bool -> unit
    = "sunml_kinsol_dls_band"

  (* Sundials < 3.0.0 *)
  external c_dls_lapack_band : 'k serial_session -> int -> int -> bool -> unit
    = "sunml_kinsol_dls_lapack_band"

  (* Sundials < 3.0.0 *)
  external c_klu
    : 'k serial_session -> 's Matrix.Sparse.sformat -> int -> int -> unit
    = "sunml_kinsol_klu_init"

  (* Sundials < 3.0.0 *)
  external c_klu_set_ordering
    : 'k serial_session -> LSI.Klu.ordering -> unit
    = "sunml_kinsol_klu_set_ordering"

  (* Sundials < 3.0.0 *)
  external c_klu_reinit
    : 'k serial_session -> int -> int -> unit
    = "sunml_kinsol_klu_reinit"

  (* Sundials < 3.0.0 *)
  external c_superlumt : 'k serial_session -> int -> int -> int -> unit
    = "sunml_kinsol_superlumt_init"

  (* Sundials < 3.0.0 *)
  external c_superlumt_set_ordering
    : 'k serial_session -> LSI.Superlumt.ordering -> unit
    = "sunml_kinsol_superlumt_set_ordering"

  (* Sundials < 3.0.0 *)
  let klu_set_ordering session ordering =
    match session.ls_callbacks with
    | SlsKluCallback _ -> c_klu_set_ordering session ordering
    | _ -> ()

  (* Sundials < 3.0.0 *)
  let klu_reinit session n onnz =
    match session.ls_callbacks with
    | SlsKluCallback _ ->
        c_klu_reinit session n (match onnz with None -> 0 | Some nnz -> nnz)
    | _ -> ()

  (* Sundials < 3.0.0 *)
  let superlumt_set_ordering session ordering =
    match session.ls_callbacks with
    | SlsSuperlumtCallback _ ->
        c_superlumt_set_ordering session ordering
    | _ -> ()

  module LSD = LSI.Direct

  (* Sundials < 3.0.0 *)
  let make_compat (type s) (type tag)
        hasjac
        (solver : (s, 'nd, 'nk, tag) LSD.solver)
        (mat : ('k, s, 'nd, 'nk) Matrix.t) session =
    match solver with
    | LSD.Dense ->
        let m, n = Matrix.(Dense.size (unwrap mat)) in
        if m <> n then raise LinearSolver.MatrixNotSquare;
        session.neqs <- n;
        c_dls_dense session hasjac
    | LSD.LapackDense ->
        let m, n = Matrix.(Dense.size (unwrap mat)) in
        if m <> n then raise LinearSolver.MatrixNotSquare;
        session.neqs <- n;
        c_dls_lapack_dense session hasjac

    | LSD.Band ->
        let open Matrix.Band in
        let { n; mu; ml } = dims (Matrix.unwrap mat) in
        session.neqs <- n;
        c_dls_band session mu ml hasjac
    | LSD.LapackBand ->
        let open Matrix.Band in
        let { n; mu; ml } = dims (Matrix.unwrap mat) in
        session.neqs <- n;
        c_dls_lapack_band session mu ml hasjac

    | LSD.Klu sinfo ->
        if not Config.klu_enabled
          then raise Config.NotImplementedBySundialsVersion;
        let smat = Matrix.unwrap mat in
        let m, n = Matrix.Sparse.size smat in
        let nnz, _ = Matrix.Sparse.dims smat in
        if m <> n then raise LinearSolver.MatrixNotSquare;
        let open LSI.Klu in
        sinfo.set_ordering <- klu_set_ordering session;
        sinfo.reinit <- klu_reinit session;
        session.neqs <- n;
        c_klu session (Matrix.Sparse.sformat smat) m nnz;
        (match sinfo.ordering with None -> ()
                                 | Some o -> c_klu_set_ordering session o)

    | LSD.Superlumt sinfo ->
        if not Config.superlumt_enabled
          then raise Config.NotImplementedBySundialsVersion;
        let smat = Matrix.unwrap mat in
        let m, n = Matrix.Sparse.size smat in
        let nnz, _ = Matrix.Sparse.dims smat in
        if m <> n then raise LinearSolver.MatrixNotSquare;
        let open LSI.Superlumt in
        sinfo.set_ordering <- superlumt_set_ordering session;
        session.neqs <- n;
        c_superlumt session m nnz sinfo.num_threads;
        (match sinfo.ordering with None -> ()
                                 | Some o -> c_superlumt_set_ordering session o)

    | LSD.Custom _ ->
        assert false

  let check_dqjac jac mat = let open Matrix in
    match get_id mat with
    | Dense | Band -> ()
    | _ -> if jac = None then invalid_arg "A Jacobian function is required"

  let set_ls_callbacks (type m) (type tag)
        ?jac (solver : (m, 'nd, 'nk, tag) LSD.solver)
        (mat : ('mk, m, 'nd, 'nk) Matrix.t) session =
    let cb = { jacfn = (match jac with None -> no_callback | Some f -> f);
               jmat  = (None : m option) } in
    begin match solver with
    | LSD.Dense ->
        session.ls_callbacks <- DlsDenseCallback cb
    | LSD.LapackDense ->
        session.ls_callbacks <- DlsDenseCallback cb
    | LSD.Band ->
        session.ls_callbacks <- DlsBandCallback cb
    | LSD.LapackBand ->
        session.ls_callbacks <- DlsBandCallback cb
    | LSD.Klu _ ->
        if jac = None then invalid_arg "Klu requires Jacobian function";
        session.ls_callbacks <- SlsKluCallback cb
    | LSD.Superlumt _ ->
        if jac = None then invalid_arg "Superlumt requires Jacobian function";
        session.ls_callbacks <- SlsSuperlumtCallback cb
    | LSD.Custom _ ->
        check_dqjac jac mat;
        session.ls_callbacks <- DirectCustomCallback cb
    end;
    session.ls_precfns <- NoPrecFns

  (* Sundials >= 3.0.0 *)
  external c_dls_set_linear_solver
    : 'k serial_session
      -> ('m, Nvector_serial.data, 'k) LSD.cptr
      -> ('mk, 'm, Nvector_serial.data, 'k) Matrix.t
      -> bool
      -> unit
    = "sunml_kinsol_dls_set_linear_solver"

  let solver ?jac ((LSD.S { LSD.rawptr; LSD.solver; LSD.matrix }) as ls)
             session nv =
    set_ls_callbacks ?jac solver matrix session;
    if in_compat_mode then make_compat (jac <> None) solver matrix session
    else c_dls_set_linear_solver session rawptr matrix (jac <> None);
    LSD.attach ls;
    session.ls_solver <- LSI.DirectSolver ls

  (* Sundials < 3.0.0 *)
  let invalidate_callback session =
    if in_compat_mode then
      match session.ls_callbacks with
      | DlsDenseCallback ({ jmat = Some d } as cb) ->
          Matrix.Dense.invalidate d;
          cb.jmat <- None
      | DlsBandCallback  ({ jmat = Some d } as cb) ->
          Matrix.Band.invalidate d;
          cb.jmat <- None
      | SlsKluCallback ({ jmat = Some d } as cb) ->
          Matrix.Sparse.invalidate d;
          cb.jmat <- None
      | SlsSuperlumtCallback ({ jmat = Some d } as cb) ->
          Matrix.Sparse.invalidate d;
          cb.jmat <- None
      | _ -> ()

  external get_work_space : 'k serial_session -> int * int
      = "sunml_kinsol_dls_get_work_space"

  let get_work_space s =
    ls_check_direct s;
    get_work_space s

  external c_get_num_jac_evals : 'k serial_session -> int
      = "sunml_kinsol_get_num_jac_evals"

  (* Sundials < 3.0.0 *)
  external c_klu_get_num_jac_evals : 'k serial_session -> int
    = "sunml_kinsol_klu_get_num_jac_evals"

  (* Sundials < 3.0.0 *)
  external c_superlumt_get_num_jac_evals : 'k serial_session -> int
    = "sunml_kinsol_superlumt_get_num_jac_evals"

  let compat_get_num_jac_evals s =
    match s.ls_callbacks with
    | SlsKluCallback _ -> c_klu_get_num_jac_evals s
    | SlsSuperlumtCallback _ -> c_superlumt_get_num_jac_evals s
    | _ -> c_get_num_jac_evals s

  let get_num_jac_evals s =
    ls_check_direct s;
    if in_compat_mode then compat_get_num_jac_evals s else
    c_get_num_jac_evals s

  external get_num_func_evals : 'k serial_session -> int
      = "sunml_kinsol_dls_get_num_func_evals"

  let get_num_func_evals s =
    ls_check_direct s;
    get_num_func_evals s
end (* }}} *)

module Spils = struct (* {{{ *)
  include SpilsTypes
  include LinearSolver.Iterative

  (* Sundials < 3.0.0 *)
  external c_spgmr : ('a, 'k) session -> int -> unit
    = "sunml_kinsol_spils_spgmr"

  (* Sundials < 3.0.0 *)
  external c_spfgmr : ('a, 'k) session -> int -> unit
    = "sunml_kinsol_spils_spfgmr"

  (* Sundials < 3.0.0 *)
  external c_spbcgs : ('a, 'k) session -> int -> unit
    = "sunml_kinsol_spils_spbcgs"

  (* Sundials < 3.0.0 *)
  external c_sptfqmr : ('a, 'k) session -> int -> unit
    = "sunml_kinsol_spils_sptfqmr"

  (* Sundials < 3.0.0 *)
  external c_set_max_restarts : ('a, 'k) session -> int -> unit
    = "sunml_kinsol_spils_set_max_restarts"

  external c_set_jac_times_vec_fn : ('a, 'k) session -> bool -> unit
    = "sunml_kinsol_spils_set_jac_times_vec_fn"

  external c_set_preconditioner
    : ('a, 'k) session -> bool -> unit
    = "sunml_kinsol_spils_set_preconditioner"

  external c_spils_set_linear_solver
    : ('a, 'k) session -> ('a, 'k) LSI.Iterative.cptr -> unit
    = "sunml_kinsol_spils_set_linear_solver"

  let old_set_max_restarts s t =
    ls_check_spils s;
    c_set_max_restarts s t

  let init_preconditioner solve setup session nv =
    c_set_preconditioner session (setup <> None);
    session.ls_precfns <- PrecFns { prec_solve_fn = solve;
                                    prec_setup_fn = setup }

  let prec_none = LSI.Iterative.(PrecNone,
                    fun session nv -> session.ls_precfns <- NoPrecFns)

  let prec_right ?setup solve = LSI.Iterative.(PrecRight,
                                              init_preconditioner solve setup)

  let not_implemented _ = raise Config.NotImplementedBySundialsVersion

  let solver (type s)
        ({ LSI.Iterative.rawptr;
           LSI.Iterative.solver;
           LSI.Iterative.compat; } as ls)
        ?jac_times_vec (prec_type, set_prec) session nv =
    if in_compat_mode then begin
      let open LSI.Iterative in
      (match (solver : ('nd, 'nk, s) solver) with
       | Spgmr ->
           c_spgmr session compat.maxl;
           (match compat.max_restarts with
            | None -> () | Some t -> c_set_max_restarts session t);
           compat.set_gs_type <- not_implemented;
           compat.set_prec_type <- not_implemented;
           compat.set_max_restarts <- old_set_max_restarts session
       | Spfgmr ->
           c_spfgmr session compat.maxl;
           (match compat.max_restarts with
            | None -> () | Some t -> c_set_max_restarts session t);
           compat.set_gs_type <- not_implemented;
           compat.set_prec_type <- not_implemented;
           compat.set_max_restarts <- old_set_max_restarts session
       | Spbcgs ->
           c_spbcgs session compat.maxl;
           compat.set_maxl <- not_implemented;
           compat.set_prec_type <- not_implemented
       | Sptfqmr ->
           c_sptfqmr session compat.maxl;
           compat.set_maxl <- not_implemented;
           compat.set_prec_type <- not_implemented
       | _ -> raise Config.NotImplementedBySundialsVersion);
      session.ls_solver <- LSI.IterativeSolver ls;
      set_prec session nv;
      session.ls_callbacks <- SpilsCallback jac_times_vec;
      if jac_times_vec <> None then c_set_jac_times_vec_fn session true
    end else
      c_spils_set_linear_solver session rawptr;
      LSI.Iterative.attach ls;
      session.ls_solver <- LSI.IterativeSolver ls;
      LSI.Iterative.(c_set_prec_type rawptr solver prec_type false);
      set_prec session nv;
      session.ls_callbacks <- SpilsCallback jac_times_vec;
      if jac_times_vec <> None then c_set_jac_times_vec_fn session true

  let set_jac_times s f =
    match s.ls_callbacks with
    | SpilsCallback _ ->
        c_set_jac_times_vec_fn s true;
        s.ls_callbacks <- SpilsCallback (Some f)
    | _ -> raise LinearSolver.InvalidLinearSolver

  let clear_jac_times s =
    match s.ls_callbacks with
    | SpilsCallback _ ->
        c_set_jac_times_vec_fn s false;
        s.ls_callbacks <- SpilsCallback None
    | _ -> raise LinearSolver.InvalidLinearSolver

  let set_preconditioner s ?setup solve =
    match s.ls_callbacks with
    | SpilsCallback _ ->
        c_set_preconditioner s (setup <> None);
        s.ls_precfns <- PrecFns { prec_setup_fn = setup;
                                  prec_solve_fn = solve }
    | _ -> raise LinearSolver.InvalidLinearSolver

  external get_work_space       : ('a, 'k) session -> int * int
      = "sunml_kinsol_spils_get_work_space"

  let get_work_space s =
    ls_check_spils s;
    get_work_space s

  external get_num_lin_iters    : ('a, 'k) session -> int
      = "sunml_kinsol_get_num_lin_iters"

  let get_num_lin_iters s =
    ls_check_spils s;
    get_num_lin_iters s

  external get_num_lin_conv_fails   : ('a, 'k) session -> int
      = "sunml_kinsol_get_num_lin_conv_fails"

  let get_num_lin_conv_fails s =
    ls_check_spils s;
    get_num_lin_conv_fails s

  external get_num_prec_evals   : ('a, 'k) session -> int
      = "sunml_kinsol_get_num_prec_evals"

  let get_num_prec_evals s =
    ls_check_spils s;
    get_num_prec_evals s

  external get_num_prec_solves  : ('a, 'k) session -> int
      = "sunml_kinsol_get_num_prec_solves"

  let get_num_prec_solves s =
    ls_check_spils s;
    get_num_prec_solves s

  external get_num_jtimes_evals : ('a, 'k) session -> int
      = "sunml_kinsol_get_num_jtimes_evals"

  let get_num_jtimes_evals s =
    ls_check_spils s;
    get_num_jtimes_evals s

  external get_num_func_evals    : ('a, 'k) session -> int
      = "sunml_kinsol_spils_get_num_func_evals"

  let get_num_func_evals s =
    ls_check_spils s;
    get_num_func_evals s
end (* }}} *)

external set_error_file : ('a, 'k) session -> Logfile.t -> unit
    = "sunml_kinsol_set_error_file"

external c_set_err_handler_fn : ('a, 'k) session -> unit
    = "sunml_kinsol_set_err_handler_fn"

let set_err_handler_fn s ferrh =
  s.errh <- ferrh;
  c_set_err_handler_fn s

external c_clear_err_handler_fn : ('a, 'k) session -> unit
    = "sunml_kinsol_clear_err_handler_fn"

let clear_err_handler_fn s =
  s.errh <- dummy_errh;
  c_clear_err_handler_fn s

external set_info_file : ('a, 'k) session -> Logfile.t -> unit
    = "sunml_kinsol_set_info_file"

external c_set_info_handler_fn : ('a, 'k) session -> unit
    = "sunml_kinsol_set_info_handler_fn"

let set_info_handler_fn s finfoh =
  s.infoh <- finfoh;
  c_set_info_handler_fn s

external c_clear_info_handler_fn : ('a, 'k) session -> unit
    = "sunml_kinsol_clear_info_handler_fn"

let clear_info_handler_fn s =
  s.infoh <- dummy_infoh;
  c_clear_info_handler_fn s

external set_print_level : ('a, 'k) session -> print_level -> unit
    = "sunml_kinsol_set_print_level"

external set_num_max_iters : ('a, 'k) session -> int -> unit
    = "sunml_kinsol_set_num_max_iters"

external set_maa : ('a, 'k) session -> int -> unit
    = "sunml_kinsol_set_maa"

external c_set_no_init_setup : ('a, 'k) session -> bool -> unit
    = "sunml_kinsol_set_no_init_setup"

let set_no_init_setup s = c_set_no_init_setup s true
let set_init_setup s = c_set_no_init_setup s false

external c_set_no_res_mon : ('a, 'k) session -> bool -> unit
    = "sunml_kinsol_set_no_res_mon"

let set_no_res_mon s = c_set_no_res_mon s true
let set_res_mon s = c_set_no_res_mon s false

external set_max_sub_setup_calls : ('a, 'k) session -> int -> unit
    = "sunml_kinsol_set_max_sub_setup_calls"

external set_max_setup_calls : ('a, 'k) session -> int -> unit
    = "sunml_kinsol_set_max_setup_calls"

external c_set_eta_form : ('a, 'k) session -> eta_choice -> unit
    = "sunml_kinsol_set_eta_form"

external c_set_eta_const_value : ('a, 'k) session -> float -> unit
    = "sunml_kinsol_set_eta_form"

external c_set_eta_params : ('a, 'k) session -> float -> float -> unit
    = "sunml_kinsol_set_eta_form"

let float_default = function None -> 0.0 | Some v -> v

let set_eta_choice s etachoice =
  match etachoice with
  | EtaChoice1 -> ()
  | EtaChoice2 { egamma; ealpha } ->
      c_set_eta_params s (float_default egamma) (float_default ealpha)
  | EtaConstant eta ->
      c_set_eta_const_value s (float_default eta);
  c_set_eta_form s etachoice

external set_res_mon_const_value : ('a, 'k) session -> float -> unit
    = "sunml_kinsol_set_res_mon_const_value"

external c_set_res_mon_params : ('a, 'k) session -> float -> float -> unit
    = "sunml_kinsol_set_res_mon_params"

let set_res_mon_params s ?omegamin ?omegamax () =
  c_set_res_mon_params s (float_default omegamin) (float_default omegamax)

external c_set_no_min_eps : ('a, 'k) session -> bool -> unit
    = "sunml_kinsol_set_no_min_eps"

let set_no_min_eps s = c_set_no_min_eps s true
let set_min_eps s = c_set_no_min_eps s false

external set_max_newton_step : ('a, 'k) session -> float -> unit
    = "sunml_kinsol_set_max_newton_step"

external set_max_beta_fails : ('a, 'k) session -> float -> unit
    = "sunml_kinsol_set_max_beta_fails"

external set_rel_err_func : ('a, 'k) session -> float -> unit
    = "sunml_kinsol_set_rel_err_func"

external set_func_norm_tol : ('a, 'k) session -> float -> unit
    = "sunml_kinsol_set_func_norm_tol"

external set_scaled_step_tol : ('a, 'k) session -> float -> unit
    = "sunml_kinsol_set_scaled_step_tol"

external c_set_constraints : ('a, 'k) session -> ('a, 'k) nvector -> unit
    = "sunml_kinsol_set_constraints"

let set_constraints s cc =
  if Sundials_configuration.safe then s.checkvec cc;
  c_set_constraints s cc

let set_sys_func s fsys =
  s.sysfn <- fsys

external get_work_space : ('a, 'k) session -> int * int
    = "sunml_kinsol_get_work_space"

external get_num_func_evals : ('a, 'k) session -> int
    = "sunml_kinsol_get_num_func_evals"

external get_num_nonlin_solv_iters : ('a, 'k) session -> int
    = "sunml_kinsol_get_num_nonlin_solv_iters"

external get_num_beta_cond_fails : ('a, 'k) session -> int
    = "sunml_kinsol_get_num_beta_cond_fails"

external get_num_backtrack_ops : ('a, 'k) session -> int
    = "sunml_kinsol_get_num_backtrack_ops"

external get_func_norm : ('a, 'k) session -> float
    = "sunml_kinsol_get_func_norm"

external get_step_length : ('a, 'k) session -> float
    = "sunml_kinsol_get_step_length"

external c_init
    : ('a, 'k) session Weak.t -> ('a, 'k) nvector -> int option -> int option
      -> (kin_mem * c_weak_ref)
    = "sunml_kinsol_init"

external c_session_finalize : ('a, 'k) session -> unit
    = "sunml_kinsol_session_finalize"

let session_finalize s =
  Dls.invalidate_callback s;
  c_session_finalize s

let init ?max_iters ?maa ?linsolv f u0 =
  let checkvec = Nvector.check u0 in
  let weakref = Weak.create 1 in
  let kin_mem, backref = c_init weakref u0 max_iters maa
  in
  let session = {
          kinsol       = kin_mem;
          backref      = backref;
          initvec      = u0;
          checkvec     = checkvec;

          exn_temp     = None;

          neqs         = 0;

          sysfn        = f;
          errh         = dummy_errh;
          infoh        = dummy_infoh;

          ls_solver    = LSI.NoSolver;
          ls_callbacks = NoCallbacks;
          ls_precfns   = NoPrecFns;
        } in
  Gc.finalise session_finalize session;
  Weak.set weakref 0 (Some session);
  (match linsolv with Some lsolver -> lsolver session u0 | None -> ());
  session

type strategy =
  | Newton            (** KIN_NONE *)
  | LineSearch        (** KIN_LINESEARCH *)
  | Picard            (** KIN_PICARD *)
  | FixedPoint        (** KIN_FP *)

type result =
  | Success           (** KIN_SUCCESS *)
  | InitialGuessOK    (** KIN_INITIAL_GUESS_OK *)
  | StoppedOnStepTol  (** KIN_STEP_LT_STPTOL *)

external c_solve : ('a, 'k) session -> ('a, 'k) nvector -> strategy
                  -> ('a, 'k) nvector -> ('a, 'k) nvector -> result
    = "sunml_kinsol_solve"

let solve s u strategy u_scale f_scale =
  if Sundials_configuration.safe then
    (s.checkvec u;
     s.checkvec u_scale;
     s.checkvec f_scale);
  if strategy <> FixedPoint && s.ls_callbacks = NoCallbacks
  then raise MissingLinearSolver
  else c_solve s u strategy u_scale f_scale

(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "sunml_kinsol_init_module"

let _ =
  c_init_module
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
      LinearSetupFailure None;
      LinearSolveFailure None;
      SystemFunctionFailure;
      FirstSystemFunctionFailure;
      RepeatedSystemFunctionFailure;
      VectorOpErr;
    |]
