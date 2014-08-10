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

include Ida_impl
type serial_session = (real_array, Nvector_serial.kind) session

(*
 * NB: The order of variant constructors and record fields is important!
 *     If these types are changed or augmented, the corresponding declarations
 *     in ida_ml.h (and code in ida_ml.c) must also be updated.
 *)

(* Solver exceptions *)
exception IllInput
exception TooClose
exception TooMuchWork
exception TooMuchAccuracy
exception ErrFailure
exception ConvergenceFailure
exception LinearInitFailure
exception LinearSetupFailure
exception LinearSolveFailure
exception ResFuncFailure
exception FirstResFuncFailure
exception RepeatedResFuncErr
exception RootFuncFailure
exception ConstraintFailure

(* Initial condition calculator exceptions *)
exception NoRecovery
exception BadEwt

(* get_dky exceptions *)
exception BadK
exception BadT
exception BadDky

let no_roots = (0, (fun _ _ _ _ -> ()))

type serial_linear_solver = (real_array, Nvector_serial.kind) linear_solver

type integrator_stats = {
    num_steps : int;
    num_res_evals : int;
    num_lin_solv_setups : int;
    num_err_test_fails : int;
    last_order : int;
    current_order : int;
    actual_init_step : float;
    last_step : float;
    current_step : float;
    current_time : float
  }

exception StopTimeReached

module VarType =
  struct
    let algebraic = 0.0
    let differential = 1.0

    type t = Algebraic | Differential
    let of_float x =
      if x = algebraic then Algebraic
      else if x = differential then Differential
      else invalid_arg ("invalid component type: " ^ string_of_float x)
    let to_float = function
      | Algebraic -> algebraic
      | Differential -> differential

    let string_of_var_type = function
      | Algebraic -> "Algebraic"
      | Differential -> "Differential"
    let string_of_float x = string_of_var_type (of_float x)
  end

let _ =
  List.iter (fun (nm, ex) -> Callback.register_exception nm ex)
  [
    ("ida_StopTimeReached",         StopTimeReached);
    ("ida_IllInput",                IllInput);
    ("ida_TooClose",                TooClose);
    ("ida_TooMuchWork",             TooMuchWork);
    ("ida_TooMuchAccuracy",         TooMuchAccuracy);
    ("ida_ErrFailure",              ErrFailure);
    ("ida_ConvergenceFailure",      ConvergenceFailure);
    ("ida_LinearInitFailure",       LinearInitFailure);
    ("ida_LinearSetupFailure",      LinearSetupFailure);
    ("ida_LinearSolveFailure",      LinearSolveFailure);
    ("ida_ResFuncFailure",          ResFuncFailure);
    ("ida_FirstResFuncFailure",     FirstResFuncFailure);
    ("ida_RepeatedResFuncErr",      RepeatedResFuncErr);
    ("ida_NoRecovery",              NoRecovery);
    ("ida_BadEwt",                  BadEwt);
    ("ida_RootFuncFailure",         RootFuncFailure);
    ("ida_ConstraintFailure",       ConstraintFailure);

    ("ida_BadK",                    BadK);
    ("ida_BadT",                    BadT);
    ("ida_BadDky",                  BadDky);
  ]

let call_resfn session t y y' res =
  let session = read_weak_ref session in
  adjust_retcode session true (session.resfn t y y') res

(* the roots function is called directly from C. *)

let call_errw session y ewt =
  let session = read_weak_ref session in
  try session.errw y ewt; 0
  with
  | Sundials.NonPositiveEwt -> -1
  | e -> (session.exn_temp <- Some e; -1)

let call_errh session details =
  let session = read_weak_ref session in
  try session.errh details
  with e ->
    prerr_endline ("Warning: error handler function raised an exception.  " ^
                   "This exception will not be propagated.")

let call_jacfn session jac j =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | DenseCallback f -> adjust_retcode session true (f jac) j
  | _ -> assert false

let call_bandjacfn session jac range j =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | BandCallback f -> adjust_retcode session true (f range jac) j
  | _ -> assert false

let call_precsolvefn session jac r z delta =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | SpilsCallback { prec_solve_fn = Some f } ->
    adjust_retcode session true (f jac r z) delta
  | _ -> assert false

let call_precsetupfn session jac =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | SpilsCallback { prec_setup_fn = Some f } ->
      adjust_retcode session true f jac
  | _ -> assert false

let call_jactimesfn session jac r z =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | SpilsCallback { jac_times_vec_fn = Some f } ->
      adjust_retcode session true (f jac r) z
  | _ -> assert false
(*
let call_linit session =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | AlternateCallback { linit = Some f } ->
      adjust_retcode session false f ()
  | _ -> assert false

let call_lsetup session convfail ypred fpred tmp =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | AlternateCallback { lsetup = Some f } ->
      adjust_retcode_and_bool session (f convfail ypred fpred) tmp
  | _ -> assert false

let call_lsolve session b weight ycur fcur =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | AlternateCallback { lsolve = f } ->
      adjust_retcode session true (f b weight ycur) fcur
  | _ -> assert false

let call_lfree session =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | AlternateCallback { lfree = Some f } -> adjust_retcode session false f ()
  | _ -> assert false
*)

let _ =
  Callback.register "c_ida_call_resfn"         call_resfn;
  Callback.register "c_ida_call_errh"          call_errh;
  Callback.register "c_ida_call_errw"          call_errw;
  Callback.register "c_ida_call_jacfn"         call_jacfn;
  Callback.register "c_ida_call_bandjacfn"     call_bandjacfn;
  Callback.register "c_ida_call_presetupfn"    call_precsetupfn;
  Callback.register "c_ida_call_presolvefn"    call_precsolvefn;
  Callback.register "c_ida_call_jactimesfn"    call_jactimesfn;

external session_finalize : ('a, 'kind) session -> unit
    = "c_ida_session_finalize"

external c_init
    : ('a, 'k) session Weak.t -> float -> ('a, 'k) nvector -> ('a, 'k) nvector
      -> (ida_mem * c_weak_ref * ida_file)
    = "c_ida_init"

external c_root_init : ('a, 'k) session -> int -> unit
    = "c_ida_root_init"

let root_init session (nroots, rootsfn) =
  c_root_init session nroots;
  session.rootsfn <- rootsfn

module Dls =
  struct
    type dense_jac_fn = Ida_impl.dense_jac_fn

    external c_dls_dense : serial_session -> int -> bool -> unit
      = "c_ida_dls_dense"

    external c_dls_lapack_dense : serial_session -> int -> bool -> unit
      = "c_ida_dls_lapack_dense"

    external c_dls_band : serial_session -> int -> int -> int -> bool -> unit
      = "c_ida_dls_band"

    external c_dls_lapack_band : serial_session -> int -> int -> int -> bool
                               -> unit
      = "c_ida_dls_lapack_band"

    external set_dense_jac_fn : serial_session -> unit
        = "c_ida_dls_set_dense_jac_fn"

    let dense jac session nv nv' =
      let neqs = Sundials.RealArray.length (Sundials.unvec nv) in
      (session.ls_callbacks <- match jac with
                               | None -> NoCallbacks
                               | Some f -> DenseCallback f);
      c_dls_dense session neqs (jac <> None)

    let lapack_dense jac session nv nv' =
      let neqs = Sundials.RealArray.length (Sundials.unvec nv) in
      (session.ls_callbacks <- match jac with
                               | None -> NoCallbacks
                               | Some f -> DenseCallback f);
      c_dls_lapack_dense session neqs (jac <> None)

    type band_jac_fn = Ida_impl.band_jac_fn

    let band p jac session nv nv' =
      let neqs = Sundials.RealArray.length (Sundials.unvec nv) in
      (session.ls_callbacks <- match jac with
                              | None -> NoCallbacks
                              | Some f -> BandCallback f);
      c_dls_band session neqs p.mupper p.mlower (jac <> None)

    let lapack_band p jac session nv nv' =
      let neqs = Sundials.RealArray.length (Sundials.unvec nv) in
      (session.ls_callbacks <- match jac with
                               | None -> NoCallbacks
                               | Some f -> BandCallback f);
      c_dls_lapack_band session neqs p.mupper p.mlower (jac <> None)

    let set_dense_jac_fn s fjacfn =
      s.ls_callbacks <- DenseCallback fjacfn;
      set_dense_jac_fn s

    external clear_dense_jac_fn : serial_session -> unit
        = "c_ida_dls_clear_dense_jac_fn"

    let clear_dense_jac_fn s =
      match s.ls_callbacks with
      | DenseCallback _ -> (s.ls_callbacks <- NoCallbacks;
                            clear_dense_jac_fn s)
      | _ -> failwith "dense linear solver not in use"

    external set_band_jac_fn : serial_session -> unit
        = "c_ida_dls_set_band_jac_fn"

    let set_band_jac_fn s fbandjacfn =
      s.ls_callbacks <- BandCallback fbandjacfn;
      set_band_jac_fn s

    external clear_band_jac_fn : serial_session -> unit
        = "c_ida_dls_clear_band_jac_fn"

    let clear_band_jac_fn s =
      match s.ls_callbacks with
      | BandCallback _ -> (s.ls_callbacks <- NoCallbacks;
                           clear_band_jac_fn s)
      | _ -> failwith "banded linear solver not in use"

    external get_work_space : serial_session -> int * int
        = "c_ida_dls_get_work_space"

    external get_num_jac_evals : serial_session -> int
        = "c_ida_dls_get_num_jac_evals"

    external get_num_res_evals : serial_session -> int
        = "c_ida_dls_get_num_res_evals"
  end

module Spils =
  struct
    type gramschmidt_type = Spils.gramschmidt_type =
      | ModifiedGS
      | ClassicalGS

    type preconditioning_type = Spils.preconditioning_type =
      | PrecNone
      | PrecLeft
      | PrecRight
      | PrecBoth

    type 'a callbacks = 'a spils_callbacks =
      {
        prec_solve_fn : (('a single_tmp, 'a) jacobian_arg -> 'a -> 'a -> float
                         -> unit) option;
        prec_setup_fn : (('a triple_tmp, 'a) jacobian_arg -> unit) option;
        jac_times_vec_fn : (('a double_tmp, 'a) jacobian_arg
                            -> 'a           (* v *)
                            -> 'a           (* Jv *)
                            -> unit) option;
      }

    let no_precond = { prec_solve_fn = None;
                       prec_setup_fn = None;
                       jac_times_vec_fn = None; }

    external c_spils_spgmr
      : ('a, 'k) session -> int -> unit
      = "c_ida_spils_spgmr"

    external c_spils_spbcg
      : ('a, 'k) session -> int -> unit
      = "c_ida_spils_spbcg"

    external c_spils_sptfqmr
      : ('a, 'k) session -> int -> unit
      = "c_ida_spils_sptfqmr"

    external c_spils_set_preconditioner
      : ('a, 'k) session -> bool -> bool -> unit
      = "c_ida_spils_set_preconditioner"

    let set_precond session cb =
      session.ls_callbacks <- SpilsCallback cb;
      c_spils_set_preconditioner session
        (cb.prec_setup_fn <> None)
        (cb.jac_times_vec_fn <> None)

    let spgmr maxl cb session nv nv' =
      let maxl = match maxl with None -> 0 | Some ml -> ml in
      c_spils_spgmr session maxl;
      set_precond session cb

    let spbcg maxl cb session nv nv' =
      let maxl = match maxl with None -> 0 | Some ml -> ml in
      c_spils_spbcg session maxl;
      set_precond session cb

    let sptfqmr maxl cb session nv nv' =
      let maxl = match maxl with None -> 0 | Some ml -> ml in
      c_spils_sptfqmr session maxl;
      set_precond session cb

    external set_preconditioner  : ('a, 'k) session -> bool -> unit
        = "c_ida_set_preconditioner"

    let set_preconditioner s fprecsetupfn fprecsolvefn =
      (match s.ls_callbacks with
       | SpilsCallback cbs ->
           s.ls_callbacks <- SpilsCallback { cbs with
                               prec_setup_fn = fprecsetupfn;
                               prec_solve_fn = Some fprecsolvefn }
       | _ -> failwith "spils solver not in use");
      set_preconditioner s (fprecsetupfn <> None)

    external set_jac_times_vec_fn : ('a, 'k) session -> unit
        = "c_ida_set_jac_times_vec_fn"

    let set_jac_times_vec_fn s fjactimesfn =
      (match s.ls_callbacks with
       | SpilsCallback cbs ->
           s.ls_callbacks <- SpilsCallback { cbs with
                               jac_times_vec_fn = Some fjactimesfn }
       | _ -> failwith "spils solver not in use");
      set_jac_times_vec_fn s

    external clear_jac_times_vec_fn : ('a, 'k) session -> unit
        = "c_ida_clear_jac_times_vec_fn"

    let clear_jac_times_vec_fn s =
      (match s.ls_callbacks with
       | SpilsCallback cbs ->
           s.ls_callbacks <- SpilsCallback { cbs with jac_times_vec_fn = None }
       | _ -> failwith "spils solver not in use");
      clear_jac_times_vec_fn s

    external set_gs_type : ('a, 'k) session -> Spils.gramschmidt_type -> unit
        = "c_ida_spils_set_gs_type"

    external set_eps_lin            : ('a, 'k) session -> float -> unit
        = "c_ida_spils_set_eps_lin"

    external set_maxl               : ('a, 'k) session -> int -> unit
        = "c_ida_spils_set_maxl"

    external get_num_lin_iters      : ('a, 'k) session -> int
        = "c_ida_spils_get_num_lin_iters"

    external get_num_conv_fails     : ('a, 'k) session -> int
        = "c_ida_spils_get_num_conv_fails"

    external get_work_space         : ('a, 'k) session -> int * int
        = "c_ida_spils_get_work_space"

    external get_num_prec_evals     : ('a, 'k) session -> int
        = "c_ida_spils_get_num_prec_evals"

    external get_num_prec_solves    : ('a, 'k) session -> int
        = "c_ida_spils_get_num_prec_solves"

    external get_num_jtimes_evals   : ('a, 'k) session -> int
        = "c_ida_spils_get_num_jtimes_evals"

    external get_num_res_evals      : ('a, 'k) session -> int
        = "c_ida_spils_get_num_res_evals"

  end
(*
module Alternate =
  struct
    type conv_fail = Cvode_impl.conv_fail =
      | NoFailures
      | FailBadJ
      | FailOther

    type 'data callbacks = 'data alternate_linsolv =
      {
        linit   : (unit -> bool) option;
        lsetup  : (conv_fail -> 'data -> 'data -> 'data triple_tmp -> bool)
                  option;
        lsolve  : 'data -> 'data -> 'data -> 'data -> unit;
        lfree   : (unit -> unit) option;
      }

    external c_set_alternate
      : ('data, 'kind) session -> bool -> bool -> bool -> unit
      = "c_ida_set_alternate"

    let make_solver f s nv =
      let { linit; lsetup; lsolve; lfree } as cb = f s nv in
      c_set_alternate s (linit <> None) (lsetup <> None) (lfree <> None);
      s.ls_callbacks <- AlternateCallback cb

  end
*)

let set_linear_solver session solver nv nv' =
  session.ls_callbacks <- NoCallbacks;
  solver session nv nv'

external sv_tolerances  : ('a, 'k) session -> float -> ('a, 'k) nvector -> unit
    = "c_ida_sv_tolerances"
external ss_tolerances  : ('a, 'k) session -> float -> float -> unit
    = "c_ida_ss_tolerances"
external wf_tolerances  : ('a, 'k) session -> unit
    = "c_ida_wf_tolerances"

type ('a, 'k) tolerance =
  | SStolerances of float * float
  | SVtolerances of float * ('a, 'k) nvector
  | WFtolerances of ('a -> 'a -> unit)

let default_tolerances = SStolerances (1.0e-4, 1.0e-8)

let set_tolerances s tol =
  match tol with
  | SStolerances (rel, abs) -> ss_tolerances s rel abs
  | SVtolerances (rel, abs) -> sv_tolerances s rel abs
  | WFtolerances ferrw -> (s.errw <- ferrw; wf_tolerances s)

let init linsolv tol resfn ?(roots=no_roots) ?(t0=0.) y y' =
  let (nroots, rootsfn) = roots in
  if nroots < 0 then
    raise (Invalid_argument "number of root functions is negative");
  (* FIXME: can we check y and y' have the same length, at least for
     some nvector types?  *)
  let weakref = Weak.create 1 in
  let ida_mem, backref, err_file = c_init weakref t0 y y' in
  (* ida_mem and backref have to be immediately captured in a session and
     associated with the finalizer before we do anything else.  *)
  let session = { ida        = ida_mem;
                  backref    = backref;
                  nroots     = nroots;
                  err_file   = err_file;
                  exn_temp   = None;
                  resfn      = resfn;
                  rootsfn    = rootsfn;
                  errh       = dummy_errh;
                  errw       = dummy_errw;
                  ls_callbacks = NoCallbacks;
                  safety_check_flags = 0;
                  sensext    = NoSensExt;
                }
  in
  Gc.finalise session_finalize session;
  Weak.set weakref 0 (Some session);
  (* Now the session is safe to use.  If any of the following fails and raises
     an exception, the GC will take care of freeing ida_mem and backref.  *)
  if nroots > 0 then
    c_root_init session nroots;
  set_linear_solver session linsolv y y';
  set_tolerances session tol;
  session

let nroots { nroots } = nroots

external c_reinit
    : ('a, 'k) session -> float -> ('a, 'k) nvector -> ('a, 'k) nvector -> unit
    = "c_ida_reinit"
let reinit session ?linsolv ?roots t0 y0 y'0 =
  c_reinit session t0 y0 y'0;
  (match linsolv with
   | None -> ()
   | Some linsolv -> set_linear_solver session linsolv y0 y'0);
  (match roots with
   | None -> ()
   | Some roots -> root_init session roots)

external get_root_info  : ('a, 'k) session -> Sundials.Roots.t -> unit
    = "c_ida_get_root_info"

external solve_normal : ('a, 'k) session -> float
                      -> ('a, 'k) nvector -> ('a,'k) nvector
                      -> float * Sundials.solver_result
    = "c_ida_solve_normal"

external solve_one_step : ('a, 'k) session -> float
                        -> ('a, 'k) nvector-> ('a, 'k) nvector
                        -> float * Sundials.solver_result
    = "c_ida_solve_one_step"

external get_dky
    : ('a, 'k) session -> float -> int -> ('a, 'k) nvector -> unit
    = "c_ida_get_dky"

external get_integrator_stats : ('a, 'k) session -> integrator_stats
    = "c_ida_get_integrator_stats"

external get_work_space         : ('a, 'k) session -> int * int
    = "c_ida_get_work_space"

external get_num_steps          : ('a, 'k) session -> int
    = "c_ida_get_num_steps"

external get_num_res_evals      : ('a, 'k) session -> int
    = "c_ida_get_num_res_evals"

external get_num_lin_solv_setups : ('a, 'k) session -> int
    = "c_ida_get_num_lin_solv_setups"

external get_num_err_test_fails : ('a, 'k) session -> int
    = "c_ida_get_num_err_test_fails"

external get_last_order         : ('a, 'k) session -> int
    = "c_ida_get_last_order"

external get_current_order      : ('a, 'k) session -> int
    = "c_ida_get_current_order"

external get_actual_init_step   : ('a, 'k) session -> float
    = "c_ida_get_actual_init_step"

external get_last_step          : ('a, 'k) session -> float
    = "c_ida_get_last_step"

external get_current_step       : ('a, 'k) session -> float
    = "c_ida_get_current_step"

external get_current_time       : ('a, 'k) session -> float
    = "c_ida_get_current_time"

let print_integrator_stats s =
  let stats = get_integrator_stats s
  in
    Printf.printf "num_steps = %d\n"           stats.num_steps;
    Printf.printf "num_res_evals = %d\n"       stats.num_res_evals;
    Printf.printf "num_lin_solv_setups = %d\n" stats.num_lin_solv_setups;
    Printf.printf "num_err_test_fails = %d\n"  stats.num_err_test_fails;
    Printf.printf "last_order = %d\n"          stats.last_order;
    Printf.printf "current_order = %d\n"       stats.current_order;
    Printf.printf "actual_init_step = %e\n"    stats.actual_init_step;
    Printf.printf "last_step = %e\n"           stats.last_step;
    Printf.printf "current_step = %e\n"        stats.current_step;
    Printf.printf "current_time = %e\n"        stats.current_time;

external set_error_file : ('a, 'k) session -> string -> bool -> unit
    = "c_ida_set_error_file"

external set_err_handler_fn  : ('a, 'k) session -> unit
    = "c_ida_set_err_handler_fn"

let set_err_handler_fn s ferrh =
  s.errh <- ferrh;
  set_err_handler_fn s

external clear_err_handler_fn  : ('a, 'k) session -> unit
    = "c_ida_clear_err_handler_fn"

let clear_err_handler_fn s =
  s.errh <- (fun _ -> ());
  clear_err_handler_fn s

external set_max_ord            : ('a, 'k) session -> int -> unit
    = "c_ida_set_max_ord"
external set_max_num_steps      : ('a, 'k) session -> int -> unit
    = "c_ida_set_max_num_steps"
external set_init_step          : ('a, 'k) session -> float -> unit
    = "c_ida_set_init_step"
external set_max_step           : ('a, 'k) session -> float -> unit
    = "c_ida_set_max_step"
external set_stop_time          : ('a, 'k) session -> float -> unit
    = "c_ida_set_stop_time"
external set_max_err_test_fails : ('a, 'k) session -> int -> unit
    = "c_ida_set_max_err_test_fails"
external set_max_nonlin_iters   : ('a, 'k) session -> int -> unit
    = "c_ida_set_max_nonlin_iters"
external set_max_conv_fails     : ('a, 'k) session -> int -> unit
    = "c_ida_set_max_conv_fails"
external set_nonlin_conv_coef   : ('a, 'k) session -> float -> unit
    = "c_ida_set_nonlin_conv_coef"

type root_direction = Sundials.RootDirs.root_direction

external set_root_direction'   : ('a, 'k) session -> Sundials.RootDirs.t -> unit
    = "c_ida_set_root_direction"

let set_root_direction s rda =
  set_root_direction' s (Sundials.RootDirs.copy_n (nroots s) rda)

let set_all_root_directions s rd =
  set_root_direction' s (Sundials.RootDirs.make (nroots s) rd)

external set_no_inactive_root_warn      : ('a, 'k) session -> unit
    = "c_ida_set_no_inactive_root_warn"
(*
   IDAGetNumStabLimOrderReds appears in the sundials 2.5.0 manual on
   p.52 but there's no such function in the implementation.  It's
   probably a leftover from earlier versions or something.

external get_num_stab_lim_order_reds    : ('a, 'k) session -> int
    = "c_ida_get_num_stab_lim_order_reds"
*)
external get_tol_scale_factor           : ('a, 'k) session -> float
    = "c_ida_get_tol_scale_factor"

external get_err_weights : ('a, 'k) session -> ('a, 'k) nvector -> unit
    = "c_ida_get_err_weights"

external get_est_local_errors : ('a, 'k) session -> ('a, 'k) nvector -> unit
    = "c_ida_get_est_local_errors"

external get_num_nonlin_solv_iters      : ('a, 'k) session -> int
    = "c_ida_get_num_nonlin_solv_iters"

external get_num_nonlin_solv_conv_fails : ('a, 'k) session -> int
    = "c_ida_get_num_nonlin_solv_conv_fails"

external get_nonlin_solv_stats          : ('a, 'k) session -> int * int
    = "c_ida_get_nonlin_solv_stats"

external get_num_g_evals                : ('a, 'k) session -> int
    = "c_ida_get_num_g_evals"

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

external set_constraints : ('a,'k) session -> ('a,'k) nvector -> unit
  = "c_ida_set_constraints"

external set_id : ('a,'k) session -> ('a,'k) nvector -> unit
  = "c_ida_set_id"

let set_var_types = set_id

external set_suppress_alg : ('a,'k) session -> bool -> unit
  = "c_ida_set_suppress_alg"

external get_num_backtrack_ops : ('a,'k) session -> int
  = "c_ida_get_num_backtrack_ops"

external c_calc_ic_y : ('a,'k) session -> ('a,'k) nvector option
                       -> float -> unit
  = "c_ida_calc_ic_y"

let calc_ic_y session ?y tout1 =
  c_calc_ic_y session y tout1

external c_calc_ic_ya_yd' :
  ('a,'k) session -> ('a,'k) nvector option -> ('a,'k) nvector option
  -> ('a,'k) nvector -> float -> unit
  = "c_ida_calc_ic_ya_ydp"

let calc_ic_ya_yd' session ?y ?y' id tout1 =
  c_calc_ic_ya_yd' session y y' id tout1
