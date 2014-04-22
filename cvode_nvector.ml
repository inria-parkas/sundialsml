(***********************************************************************)
(*                                                                     *)
(*     OCaml interface to Sundials (serial) CVODE and IDA solvers      *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2013 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

include Cvode

type 'a nvector = 'a Nvector.nvector

type root_array = Sundials.Roots.t
type root_val_array = Sundials.Roots.val_array


type 'a single_tmp = 'a
type 'a triple_tmp = 'a * 'a * 'a

type ('t, 'a) jacobian_arg =
  {
    jac_t   : float;
    jac_y   : 'a;
    jac_fy  : 'a;
    jac_tmp : 't
  }

type 'a prec_solve_arg =
  {
    rhs   : 'a;
    gamma : float;
    delta : float;
    left  : bool;
  }

(* Note this definition differs from the one in cvode_serial, so the tag
   values are different.  *)
type 'a linear_solver =
  | Diag
  | Spgmr of spils_params * 'a spils_callbacks
  | Spbcg of spils_params * 'a spils_callbacks
  | Sptfqmr of spils_params * 'a spils_callbacks
  | BandedSpgmr of spils_params * bandrange
  | BandedSpbcg of spils_params * bandrange
  | BandedSptfqmr of spils_params * bandrange
and bandrange = { mupper : int;
                  mlower : int; }
and spils_params = { maxl : int;
                     prec_type : preconditioning_type; }
and 'a spils_callbacks =
  {
    prec_solve_fn : (('a single_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                     -> 'a nvector -> unit) option;
    prec_setup_fn : (('a triple_tmp, 'a) jacobian_arg -> bool -> float -> bool)
                    option;
    jac_times_vec_fn :
      (('a single_tmp, 'a) jacobian_arg
       -> 'a (* v *)
       -> 'a (* Jv *)
       -> unit) option;
  }
let spils_no_precond = { prec_solve_fn = None;
                         prec_setup_fn = None;
                         jac_times_vec_fn = None; }

type 'a iter =
  | Newton of 'a linear_solver
  | Functional

type cvode_mem
type cvode_file
type c_weak_ref
type 'a session = {
        cvode      : cvode_mem;
        backref    : c_weak_ref;
        neqs       : int;
        nroots     : int;
        err_file   : cvode_file;

        mutable exn_temp   : exn option;

        mutable rhsfn      : float -> 'a -> 'a -> unit;
        mutable rootsfn    : float -> 'a -> root_val_array -> unit;
        mutable errh       : error_details -> unit;
        mutable errw       : 'a -> 'a -> unit;
        mutable jacfn      : ('a triple_tmp, 'a) jacobian_arg -> Dls.DenseMatrix.t -> unit;
        mutable bandjacfn  : ('a triple_tmp, 'a) jacobian_arg -> int -> int
                               -> Dls.BandMatrix.t -> unit;
        mutable presetupfn : ('a triple_tmp, 'a) jacobian_arg -> bool -> float -> bool;
        mutable presolvefn : ('a single_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                               -> 'a nvector -> unit;
        mutable jactimesfn : ('a single_tmp, 'a) jacobian_arg -> 'a -> 'a -> unit;
      }

(* interface *)

let read_weak_ref x : 'a session =
  match Weak.get x 0 with
  | Some y -> y
  | None -> raise (Failure "Internal error: weak reference is dead")
let adjust_retcode = fun session check_recoverable f x ->
  try f x; 0
  with
  | RecoverableFailure when check_recoverable -> 1
  | e -> (session.exn_temp <- Some e; -1)
let call_rhsfn session t y y' =
  let session = read_weak_ref session in
  adjust_retcode session true (session.rhsfn t y) y'
let call_errw session y ewt =
  let session = read_weak_ref session in
  adjust_retcode session false (session.errw y) ewt
let call_errh session details =
  let session = read_weak_ref session in
  try session.errh details
  with e ->
    prerr_endline ("Warning: error handler function raised an exception.  " ^
                   "This exception will not be propagated.")
let call_jacfn session jac j =
  let session = read_weak_ref session in
  adjust_retcode session true (session.jacfn jac) j
let call_bandjacfn session jac mupper mlower j =
  let session = read_weak_ref session in
  adjust_retcode session true (session.bandjacfn jac mupper mlower) j
let call_presolvefn session jac r z =
  let session = read_weak_ref session in
  adjust_retcode session true (session.presolvefn jac r) z
let call_jactimesfn session jac v jv =
  let session = read_weak_ref session in
  adjust_retcode session true (session.jactimesfn jac v) jv
let _ =
  Callback.register "c_nvec_cvode_call_rhsfn"         call_rhsfn;
  Callback.register "c_nvec_cvode_call_errh"          call_errh;
  Callback.register "c_nvec_cvode_call_errw"          call_errw;
  Callback.register "c_nvec_cvode_call_jacfn"         call_jacfn;
  Callback.register "c_nvec_cvode_call_bandjacfn"     call_bandjacfn;
  Callback.register "c_nvec_cvode_call_presolvefn"    call_presolvefn;
  Callback.register "c_nvec_cvode_call_jactimesfn"    call_jactimesfn;

external session_finalize : 'a session -> unit
    = "c_cvode_session_finalize"

external c_init
    : 'a session Weak.t -> lmm -> 'a iter -> 'a nvector -> float
      -> (cvode_mem * c_weak_ref * cvode_file)
    = "c_nvec_cvode_init"

external c_root_init : 'a session -> int -> unit
    = "c_nvec_cvode_root_init"
let root_init session (nroots, rootsfn) =
  c_root_init session nroots;
  session.rootsfn <- rootsfn

external c_diag : 'a session -> unit
  = "c_cvode_diag"

external c_spils_spgmr
  : 'a session -> int -> preconditioning_type -> unit
  = "c_nvec_cvode_spils_spgmr"

external c_spils_spbcg
  : 'a session -> int -> preconditioning_type -> unit
  = "c_nvec_cvode_spils_spbcg"

external c_spils_sptfqmr
  : 'a session -> int -> preconditioning_type -> unit
  = "c_nvec_cvode_spils_sptfqmr"

external c_spils_banded_spgmr
  : 'a session -> int -> int -> int -> preconditioning_type -> unit
  = "c_nvec_cvode_spils_banded_spgmr"

external c_spils_banded_spbcg
  : 'a session -> int -> int -> int -> preconditioning_type -> unit
  = "c_nvec_cvode_spils_banded_spbcg"

external c_spils_banded_sptfqmr
  : 'a session -> int -> int -> int -> preconditioning_type -> unit
  = "c_nvec_cvode_spils_banded_sptfqmr"

external c_spils_set_preconditioner
  : 'a session -> bool -> bool -> unit
  = "c_ba_cvode_spils_set_preconditioner"

external c_set_functional : 'a session -> unit
  = "c_cvode_set_functional"

let shouldn't_be_called fcn =
  failwith ("internal error in sundials: " ^ fcn ^ " is called")
let dummy_dense_jac _ _ = shouldn't_be_called "dummy_dense_jac"
let dummy_band_jac _ _ _ _ = shouldn't_be_called "dummy_band_jac"
let dummy_prec_setup _ _ _ = shouldn't_be_called "dummy_prec_setup"
let dummy_prec_solve _ _ _ = shouldn't_be_called "dummy_prec_solve"
let dummy_jac_times_vec _ _ _ = shouldn't_be_called "dummy_jac_times_vec"

let set_iter_type session iter =
  let optionally f = function
    | None -> ()
    | Some x -> f x
  in
  let set_precond prec_type cb =
    match prec_type with
    | PrecNone -> ()
    | PrecLeft | PrecRight | PrecBoth ->
      match cb.prec_solve_fn with
      | None -> invalid_arg "preconditioning type is not PrecNone, but no \
                             solve function given"
      | Some solve_fn ->
        c_spils_set_preconditioner session
          (cb.prec_setup_fn <> None)
          (cb.jac_times_vec_fn <> None);
        session.presolvefn <- solve_fn;
        optionally (fun f -> session.presetupfn <- f) cb.prec_setup_fn;
        optionally (fun f -> session.jactimesfn <- f) cb.jac_times_vec_fn
  in
  (* Release references to all linear solver-related callbacks.  *)
  session.jacfn      <- dummy_dense_jac;
  session.bandjacfn  <- dummy_band_jac;
  session.presetupfn <- dummy_prec_setup;
  session.presolvefn <- dummy_prec_solve;
  session.jactimesfn <- dummy_jac_times_vec;
  match iter with
  | Functional -> c_set_functional session
  | Newton linsolv ->
    (* Iter type will be set to CV_NEWTON in the functions that set the linear
       solver.  *)
    match linsolv with
    | Diag -> c_diag session
    | Spgmr (par, cb) ->
      c_spils_spgmr session par.maxl par.prec_type;
      set_precond par.prec_type cb
    | Spbcg (par, cb) ->
      c_spils_spbcg session par.maxl par.prec_type;
      set_precond par.prec_type cb
    | Sptfqmr (par, cb) ->
      c_spils_sptfqmr session par.maxl par.prec_type;
      set_precond par.prec_type cb
    | BandedSpgmr (sp, br) ->
      c_spils_banded_spgmr session br.mupper br.mlower sp.maxl sp.prec_type
    | BandedSpbcg (sp, br) ->
      c_spils_banded_spbcg session br.mupper br.mlower sp.maxl sp.prec_type
    | BandedSptfqmr (sp, br) ->
      c_spils_banded_sptfqmr session br.mupper br.mlower sp.maxl sp.prec_type

external sv_tolerances  : 'a session -> float -> 'a nvector -> unit
    = "c_nvec_cvode_sv_tolerances"
external ss_tolerances  : 'a session -> float -> float -> unit
    = "c_cvode_ss_tolerances"
external wf_tolerances  : 'a session -> unit
    = "c_nvec_cvode_wf_tolerances"

let init lmm iter f ?(roots=no_roots) ?(t0=0.) (neqs, y0) =
  let (nroots, roots) = roots in
  let weakref = Weak.create 1 in
  let cvode_mem, backref, err_file = c_init weakref lmm iter y0 t0 in
  (* cvode_mem and backref have to be immediately captured in a session and
     associated with the finalizer before we do anything else.  *)
  let session = {
          cvode      = cvode_mem;
          backref    = backref;
          neqs       = neqs;
          nroots     = nroots;
          err_file   = err_file;

          exn_temp   = None;

          rhsfn      = f;
          rootsfn    = roots;
          errh       = (fun _ -> ());
          errw       = (fun _ _ -> ());
          jacfn      = dummy_dense_jac;
          bandjacfn  = dummy_band_jac;
          presetupfn = dummy_prec_setup;
          presolvefn = dummy_prec_solve;
          jactimesfn = dummy_jac_times_vec;
        } in
  Gc.finalise session_finalize session;
  Weak.set weakref 0 (Some session);
  (* Now the session is safe to use.  If any of the following fails and raises
     an exception, the GC will take care of freeing cvode_mem and backref.  *)
  if nroots > 0 then
    c_root_init session nroots;
  set_iter_type session iter;
  ss_tolerances session 1.0e-4 1.0e-8;
  session

let nroots { nroots } = nroots
let neqs { neqs } = neqs

external c_reinit
    : 'a session -> float -> 'a nvector -> unit
    = "c_nvec_cvode_reinit"
let reinit session ?iter_type ?roots t0 y0 =
  c_reinit session t0 y0;
  (match iter_type with
   | None -> ()
   | Some iter_type -> set_iter_type session iter_type);
  (match roots with
   | None -> ()
   | Some roots -> root_init session roots)

let wf_tolerances s ferrw =
  s.errw <- ferrw;
  wf_tolerances s

external get_root_info  : 'a session -> Roots.t -> unit
    = "c_cvode_get_root_info"

external solve_normal
    : 'a session -> float -> 'a nvector -> float * solver_result
    = "c_nvec_cvode_solve_normal"

external solve_one_step
    : 'a session -> float -> 'a nvector -> float * solver_result
    = "c_nvec_cvode_solve_one_step"

external get_dky
    : 'a session -> float -> int -> 'a nvector -> unit
    = "c_nvec_cvode_get_dky"

external get_integrator_stats : 'a session -> integrator_stats
    = "c_cvode_get_integrator_stats"

external get_work_space         : 'a session -> int * int
    = "c_cvode_get_work_space"

external get_num_steps          : 'a session -> int
    = "c_cvode_get_num_steps"

external get_num_rhs_evals      : 'a session -> int
    = "c_cvode_get_num_rhs_evals"

external get_num_lin_solv_setups : 'a session -> int
    = "c_cvode_get_num_lin_solv_setups"

external get_num_err_test_fails : 'a session -> int
    = "c_cvode_get_num_err_test_fails"

external get_last_order         : 'a session -> int
    = "c_cvode_get_last_order"

external get_current_order      : 'a session -> int
    = "c_cvode_get_current_order"

external get_actual_init_step   : 'a session -> float
    = "c_cvode_get_actual_init_step"

external get_last_step          : 'a session -> float
    = "c_cvode_get_last_step"

external get_current_step       : 'a session -> float
    = "c_cvode_get_current_step"

external get_current_time       : 'a session -> float
    = "c_cvode_get_current_time"

let print_integrator_stats s =
  let stats = get_integrator_stats s
  in
    Printf.printf "num_steps = %d\n"           stats.num_steps;
    Printf.printf "num_rhs_evals = %d\n"       stats.num_rhs_evals;
    Printf.printf "num_lin_solv_setups = %d\n" stats.num_lin_solv_setups;
    Printf.printf "num_err_test_fails = %d\n"  stats.num_err_test_fails;
    Printf.printf "last_order = %d\n"          stats.last_order;
    Printf.printf "current_order = %d\n"       stats.current_order;
    Printf.printf "actual_init_step = %e\n"    stats.actual_init_step;
    Printf.printf "last_step = %e\n"           stats.last_step;
    Printf.printf "current_step = %e\n"        stats.current_step;
    Printf.printf "current_time = %e\n"        stats.current_time;

external set_error_file : 'a session -> string -> bool -> unit
    = "c_cvode_set_error_file"

external set_err_handler_fn  : 'a session -> unit
    = "c_nvec_cvode_set_err_handler_fn"

let set_err_handler_fn s ferrh =
  s.errh <- ferrh;
  set_err_handler_fn s

external clear_err_handler_fn  : 'a session -> unit
    = "c_nvec_cvode_clear_err_handler_fn"

let clear_err_handler_fn s =
  s.errh <- (fun _ -> ());
  clear_err_handler_fn s

external set_max_ord            : 'a session -> int -> unit
    = "c_cvode_set_max_ord"
external set_max_num_steps      : 'a session -> int -> unit
    = "c_cvode_set_max_num_steps"
external set_max_hnil_warns     : 'a session -> int -> unit
    = "c_cvode_set_max_hnil_warns"
external set_stab_lim_det       : 'a session -> bool -> unit
    = "c_cvode_set_stab_lim_det"
external set_init_step          : 'a session -> float -> unit
    = "c_cvode_set_init_step"
external set_min_step           : 'a session -> float -> unit
    = "c_cvode_set_min_step"
external set_max_step           : 'a session -> float -> unit
    = "c_cvode_set_max_step"
external set_stop_time          : 'a session -> float -> unit
    = "c_cvode_set_stop_time"
external set_max_err_test_fails : 'a session -> int -> unit
    = "c_cvode_set_max_err_test_fails"
external set_max_nonlin_iters   : 'a session -> int -> unit
    = "c_cvode_set_max_nonlin_iters"
external set_max_conv_fails     : 'a session -> int -> unit
    = "c_cvode_set_max_conv_fails"
external set_nonlin_conv_coef   : 'a session -> float -> unit
    = "c_cvode_set_nonlin_conv_coef"

external set_root_direction'    : 'a session -> RootDirs.t -> unit
    = "c_cvode_set_root_direction"

let set_root_direction s rda =
  set_root_direction' s (RootDirs.copy_n (nroots s) rda)

let set_all_root_directions s rd =
  set_root_direction' s (RootDirs.make (nroots s) rd)

external set_no_inactive_root_warn      : 'a session -> unit
    = "c_cvode_set_no_inactive_root_warn"

external get_num_stab_lim_order_reds    : 'a session -> int
    = "c_cvode_get_num_stab_lim_order_reds"

external get_tol_scale_factor           : 'a session -> float
    = "c_cvode_get_tol_scale_factor"

external get_err_weights                : 'a session -> 'a nvector -> unit
    = "c_nvec_cvode_get_err_weights"

external get_est_local_errors           : 'a session -> 'a nvector -> unit
    = "c_nvec_cvode_get_est_local_errors"

external get_num_nonlin_solv_iters      : 'a session -> int
    = "c_cvode_get_num_nonlin_solv_iters"

external get_num_nonlin_solv_conv_fails : 'a session -> int
    = "c_cvode_get_num_nonlin_solv_conv_fails"

external get_num_g_evals                : 'a session -> int
    = "c_cvode_get_num_g_evals"

module Dls =
  struct
    external set_dense_jac_fn  : 'a session -> unit
        = "c_ba_cvode_dls_set_dense_jac_fn"

    let set_dense_jac_fn s fjacfn =
      s.jacfn <- fjacfn;
      set_dense_jac_fn s

    external clear_dense_jac_fn : 'a session -> unit
        = "c_ba_cvode_dls_clear_dense_jac_fn"

    let clear_dense_jac_fn s =
      s.jacfn <- dummy_dense_jac;
      clear_dense_jac_fn s

    external set_band_jac_fn   : 'a session -> unit
        = "c_ba_cvode_dls_set_band_jac_fn"

    let set_band_jac_fn s fbandjacfn =
      s.bandjacfn <- fbandjacfn;
      set_band_jac_fn s

    external clear_band_jac_fn : 'a session -> unit
        = "c_ba_cvode_dls_clear_band_jac_fn"

    let clear_band_jac_fn s =
      s.bandjacfn <- dummy_band_jac;
      clear_band_jac_fn s

    external get_work_space : 'a session -> int * int
        = "c_cvode_dls_get_work_space"

    external get_num_jac_evals    : 'a session -> int
        = "c_cvode_dls_get_num_jac_evals"

    external get_num_rhs_evals    : 'a session -> int
        = "c_cvode_dls_get_num_rhs_evals"
  end

module Diag =
  struct
    external get_work_space       : 'a session -> int * int
        = "c_cvode_diag_get_work_space"

    external get_num_rhs_evals    : 'a session -> int
        = "c_cvode_diag_get_num_rhs_evals"
  end

module BandPrec =
  struct
    external get_work_space : 'a session -> int * int
        = "c_cvode_bandprec_get_work_space"

    external get_num_rhs_evals    : 'a session -> int
        = "c_cvode_bandprec_get_num_rhs_evals"
  end

module Spils =
  struct
    type 'a solve_arg = 'a prec_solve_arg =
      {
        rhs : 'a;
        gamma : float;
        delta : float;
        left : bool;
      }

    type gramschmidt_type =
      | ModifiedGS
      | ClassicalGS

    external set_preconditioner  : 'a session -> unit
        = "c_nvec_cvode_set_preconditioner"

    let set_preconditioner s fpresetupfn fpresolvefn =
      s.presetupfn <- fpresetupfn;
      s.presolvefn <- fpresolvefn;
      set_preconditioner s

    external set_jac_times_vec_fn : 'a session -> unit
        = "c_nvec_cvode_set_jac_times_vec_fn"

    let set_jac_times_vec_fn s fjactimesfn =
      s.jactimesfn <- fjactimesfn;
      set_jac_times_vec_fn s

    external clear_jac_times_vec_fn : 'a session -> unit
        = "c_nvec_cvode_clear_jac_times_vec_fn"

    let clear_jac_times_vec_fn s =
      s.jactimesfn <- (fun _ _ _ -> ());
      clear_jac_times_vec_fn s

    external set_prec_type : 'a session -> preconditioning_type -> unit
        = "c_cvode_set_prec_type"

    external set_gs_type : 'a session -> gramschmidt_type -> unit
        = "c_cvode_set_gs_type"

    external set_eps_lin            : 'a session -> float -> unit
        = "c_cvode_set_eps_lin"

    external set_maxl               : 'a session -> int -> unit
        = "c_cvode_set_maxl"

    external get_num_lin_iters      : 'a session -> int
        = "c_cvode_spils_get_num_lin_iters"

    external get_num_conv_fails     : 'a session -> int
        = "c_cvode_spils_get_num_conv_fails"

    external get_work_space         : 'a session -> int * int
        = "c_cvode_spils_get_work_space"

    external get_num_prec_evals     : 'a session -> int
        = "c_cvode_spils_get_num_prec_evals"

    external get_num_prec_solves    : 'a session -> int
        = "c_cvode_spils_get_num_prec_solves"

    external get_num_jtimes_evals   : 'a session -> int
        = "c_cvode_spils_get_num_jtimes_evals"

    external get_num_rhs_evals      : 'a session -> int
        = "c_cvode_spils_get_num_rhs_evals"

  end

