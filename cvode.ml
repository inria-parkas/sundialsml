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
include Cvode_impl

(*
 * NB: The order of variant constructors and record fields is important!
 *     If these types are changed or augmented, the corresponding declarations
 *     in cvode_serial.h (and code in cvode_serial.c) must also be updated.
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
exception RhsFuncFailure
exception FirstRhsFuncFailure
exception RepeatedRhsFuncFailure
exception UnrecoverableRhsFuncFailure
exception RootFuncFailure

(* get_dky exceptions *)
exception BadK
exception BadT
exception BadDky

let no_roots = (0, dummy_rootsfn)

type lmm =
  | Adams
  | BDF

type ('a, 'kind) iter =
  | Newton of ('a, 'kind) linear_solver
  | Functional

type integrator_stats = {
    num_steps : int;
    num_rhs_evals : int;
    num_lin_solv_setups : int;
    num_err_test_fails : int;
    last_order : int;
    current_order : int;
    actual_init_step : float;
    last_step : float;
    current_step : float;
    current_time : float
  }

external c_root_init : ('a, 'k) session -> int -> unit
    = "c_cvode_root_init"

let root_init session (nroots, rootsfn) =
  c_root_init session nroots;
  session.rootsfn <- rootsfn

module Diag =
  struct
    external c_cvode_diag : ('a, 'k) session -> unit
      = "c_cvode_diag"

    let solver session nv = c_cvode_diag session

    external get_work_space       : ('a, 'k) session -> int * int
        = "c_cvode_diag_get_work_space"

    external get_num_rhs_evals    : ('a, 'k) session -> int
        = "c_cvode_diag_get_num_rhs_evals"
  end

module Dls =
  struct
    include DlsTypes

    external c_dls_dense : serial_session -> int -> bool -> unit
      = "c_cvode_dls_dense"

    external c_dls_lapack_dense : serial_session -> int -> bool -> unit
      = "c_cvode_dls_lapack_dense"

    external c_dls_band : (serial_session * int) -> int -> int -> bool -> unit
      = "c_cvode_dls_band"

    external c_dls_lapack_band
      : (serial_session * int) -> int -> int -> bool -> unit
      = "c_cvode_dls_lapack_band"

    external set_dense_jac_fn : serial_session -> unit
        = "c_cvode_dls_set_dense_jac_fn"

    let dense jac session nv =
      let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
      (session.ls_callbacks <- match jac with
                               | None -> NoCallbacks
                               | Some f -> DenseCallback
                                             { jacfn = f; dmat = None });
      c_dls_dense session neqs (jac <> None)

    let lapack_dense jac session nv =
      let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
      (session.ls_callbacks <- match jac with
                               | None -> NoCallbacks
                               | Some f -> DenseCallback
                                             { jacfn = f; dmat = None });
      c_dls_lapack_dense session neqs (jac <> None)

    let band p jac session nv =
      let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
      (session.ls_callbacks <- match jac with
                               | None -> NoCallbacks
                               | Some f -> BandCallback
                                             { bjacfn = f; bmat = None });
      c_dls_band (session, neqs) p.mupper p.mlower (jac <> None)

    let lapack_band p jac session nv =
      let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
      (session.ls_callbacks <- match jac with
                               | None -> NoCallbacks
                               | Some f -> BandCallback
                                             { bjacfn = f; bmat = None });
      c_dls_lapack_band (session, neqs) p.mupper p.mlower (jac <> None)

    let invalidate_callback (type d) (type k) (session : (d, k) session) =
      match session.ls_callbacks with
      | DenseCallback ({ dmat = Some d } as cb) ->
          Dls.DenseMatrix.invalidate d;
          cb.dmat <- None
      | BandCallback  ({ bmat = Some d } as cb) ->
          Dls.BandMatrix.invalidate d;
          cb.bmat <- None
      | _ -> ()

    let set_dense_jac_fn s fjacfn =
      invalidate_callback s;
      s.ls_callbacks <- DenseCallback { jacfn = fjacfn; dmat = None };
      set_dense_jac_fn s

    external clear_dense_jac_fn : serial_session -> unit
        = "c_cvode_dls_clear_dense_jac_fn"

    let clear_dense_jac_fn s =
      match s.ls_callbacks with
      | DenseCallback _ -> (invalidate_callback s;
                            s.ls_callbacks <- NoCallbacks;
                            clear_dense_jac_fn s)
      | _ -> failwith "dense linear solver not in use"

    external set_band_jac_fn : serial_session -> unit
        = "c_cvode_dls_set_band_jac_fn"

    let set_band_jac_fn s f =
      invalidate_callback s;
      s.ls_callbacks <- BandCallback { bjacfn = f; bmat = None };
      set_band_jac_fn s

    external clear_band_jac_fn : serial_session -> unit
        = "c_cvode_dls_clear_band_jac_fn"

    let clear_band_jac_fn s =
      match s.ls_callbacks with
      | BandCallback _ -> (invalidate_callback s;
                           s.ls_callbacks <- NoCallbacks;
                           clear_band_jac_fn s)
      | _ -> failwith "banded linear solver not in use"

    external get_work_space : serial_session -> int * int
        = "c_cvode_dls_get_work_space"

    external get_num_jac_evals : serial_session -> int
        = "c_cvode_dls_get_num_jac_evals"

    external get_num_rhs_evals : serial_session -> int
        = "c_cvode_dls_get_num_rhs_evals"
  end

module Spils =
  struct
    include SpilsTypes

    external c_spgmr
      : ('a, 'k) session -> int -> Spils.preconditioning_type -> unit
      = "c_cvode_spils_spgmr"

    external c_spbcg
      : ('a, 'k) session -> int -> Spils.preconditioning_type -> unit
      = "c_cvode_spils_spbcg"

    external c_sptfqmr
      : ('a, 'k) session -> int -> Spils.preconditioning_type -> unit
      = "c_cvode_spils_sptfqmr"

    external c_set_preconditioner
      : ('a, 'k) session -> bool -> unit
      = "c_cvode_spils_set_preconditioner"

    external c_set_jac_times_vec_fn
      : ('a, 'k) session -> bool -> unit
      = "c_cvode_spils_set_jac_times_vec_fn"

    let init_preconditioner solve setup jac_times session nv =
      c_set_preconditioner session (setup <> None);
      c_set_jac_times_vec_fn session (jac_times <> None);
      session.ls_callbacks <-
        SpilsCallback { prec_solve_fn = solve;
                        prec_setup_fn = setup;
                        jac_times_vec_fn = jac_times }

    let prec_none = InternalPrecNone
    let prec_left ?setup ?jac_times_vec solve =
      InternalPrecLeft (init_preconditioner solve setup jac_times_vec)
    let prec_right ?setup ?jac_times_vec solve =
      InternalPrecRight (init_preconditioner solve setup jac_times_vec)
    let prec_both ?setup ?jac_times_vec solve =
      InternalPrecBoth (init_preconditioner solve setup jac_times_vec)

    let init_spils init maxl prec session nv =
      let with_prec prec_type set_prec =
        init session maxl prec_type;
        set_prec session nv
      in
      match prec with
      | InternalPrecNone -> init session maxl Spils.PrecNone
      | InternalPrecLeft set_prec  -> with_prec Spils.PrecLeft set_prec
      | InternalPrecRight set_prec -> with_prec Spils.PrecRight set_prec
      | InternalPrecBoth set_prec  -> with_prec Spils.PrecBoth set_prec

    let spgmr ?(maxl=0) prec session nv =
      init_spils c_spgmr maxl prec session nv

    let spbcg ?(maxl=0) prec session nv =
      init_spils c_spbcg maxl prec session nv

    let sptfqmr ?(maxl=0) prec session nv =
      init_spils c_sptfqmr maxl prec session nv

    let set_preconditioner s ?setup solve =
      match s.ls_callbacks with
      | SpilsCallback cbs ->
        c_set_preconditioner s (setup <> None);
        s.ls_callbacks <- SpilsCallback { cbs with
                                          prec_setup_fn = setup;
                                          prec_solve_fn = solve }
      | SpilsBandCallback _ -> failwith "User-defined preconditioner not in use"
      | _ -> failwith "spils solver not in use"

    let set_jac_times_vec_fn s f =
      match s.ls_callbacks with
      | SpilsCallback cbs ->
        c_set_jac_times_vec_fn s true;
        s.ls_callbacks <- SpilsCallback { cbs with jac_times_vec_fn = Some f }
      | SpilsBandCallback _ -> failwith "User-defined preconditioner not in use"
      | _ -> failwith "spils solver not in use"

    let clear_jac_times_vec_fn s =
      match s.ls_callbacks with
      | SpilsCallback cbs ->
        c_set_jac_times_vec_fn s false;
        s.ls_callbacks <- SpilsCallback { cbs with jac_times_vec_fn = None }
      | SpilsBandCallback _ -> failwith "User-defined preconditioner not in use"
      | _ -> failwith "spils solver not in use"

    external c_set_prec_type
        : ('a, 'k) session -> Spils.preconditioning_type -> unit
        = "c_cvode_spils_set_prec_type"

    let set_prec_type s = function
      | PrecNone -> c_set_prec_type s Spils.PrecNone
      | PrecLeft -> c_set_prec_type s Spils.PrecLeft
      | PrecRight -> c_set_prec_type s Spils.PrecRight
      | PrecBoth -> c_set_prec_type s Spils.PrecBoth

    external set_gs_type : ('a, 'k) session -> gramschmidt_type -> unit
        = "c_cvode_spils_set_gs_type"

    external set_eps_lin            : ('a, 'k) session -> float -> unit
        = "c_cvode_spils_set_eps_lin"

    external set_maxl               : ('a, 'k) session -> int -> unit
        = "c_cvode_spils_set_maxl"

    external get_num_lin_iters      : ('a, 'k) session -> int
        = "c_cvode_spils_get_num_lin_iters"

    external get_num_conv_fails     : ('a, 'k) session -> int
        = "c_cvode_spils_get_num_conv_fails"

    external get_work_space         : ('a, 'k) session -> int * int
        = "c_cvode_spils_get_work_space"

    external get_num_prec_evals     : ('a, 'k) session -> int
        = "c_cvode_spils_get_num_prec_evals"

    external get_num_prec_solves    : ('a, 'k) session -> int
        = "c_cvode_spils_get_num_prec_solves"

    external get_num_jtimes_evals   : ('a, 'k) session -> int
        = "c_cvode_spils_get_num_jtimes_evals"

    external get_num_rhs_evals      : ('a, 'k) session -> int
        = "c_cvode_spils_get_num_rhs_evals"

    module Banded = struct
      external c_set_preconditioner
        : ('a, 'k) session -> int -> int -> int -> unit
        = "c_cvode_spils_set_banded_preconditioner"
      (* Note: CVBandPrecInit seems to be designed only to be called on
         a fresh spils solver (i.e. right after CVSpgmr, CVSpbcg, or
         CVSptfqmr).

         As of Sundials 2.5.0, calling

           CVSpgmr -> CVBandPrecInit -> CVBandPrecInit

         triggers a memory leak.  Calling CVodeReInit does NOT help.
         The only way to prevent leakage is to allocate a fresh spils
         instance, thus:

           CVSpgmr -> CVBandPrecInit -> CVSpgmr -> CVBandPrecInit.

         If you call

           CVSpgmr -> CVSpilsSetPreconditioner -> CVBandPrecInit,

         nothing grave happens, but the memory associated with
         CVBandPrecInit won't be freed until the spils solver is torn
         down.  If you call BandPrecInit -> SetPreconditioner ->
         BandPrecInit, you also get a memory leak.

         (Perhaps set_preconditioner should be hidden too?  In that
          case, we should somehow strip set_prec_type of the ability
          to change PrecNone to anything else.)

         set_jac_times_vec_fn, as well as similar functions for Dls
         solvers, are kept because they accept NULL to remove the
         callback.  This design clearly anticipates being called
         multiple times on the same solver instance.  *)

      let init_preconditioner jac_times_vec bandrange session nv =
        c_set_preconditioner session (RealArray.length (Nvector.unwrap nv))
          bandrange.mupper bandrange.mlower;
        c_set_jac_times_vec_fn session (jac_times_vec <> None);
        session.ls_callbacks <- SpilsBandCallback jac_times_vec

      let prec_none = InternalPrecNone
      let prec_left ?jac_times_vec bandrange =
        InternalPrecLeft (init_preconditioner jac_times_vec bandrange)
      let prec_right ?jac_times_vec bandrange =
        InternalPrecRight (init_preconditioner jac_times_vec bandrange)
      let prec_both ?jac_times_vec bandrange =
        InternalPrecBoth (init_preconditioner jac_times_vec bandrange)

      let set_jac_times_vec_fn s f =
        match s.ls_callbacks with
        | SpilsCallback _ -> failwith "Banded preconditioner not in use"
        | SpilsBandCallback _ ->
          c_set_jac_times_vec_fn s true;
          s.ls_callbacks <- SpilsBandCallback (Some f)
        | _ -> failwith "spils solver not in use"

      let clear_jac_times_vec_fn s =
        match s.ls_callbacks with
        | SpilsCallback _ -> failwith "Banded preconditioner not in use"
        | SpilsBandCallback _ ->
          c_set_jac_times_vec_fn s false;
          s.ls_callbacks <- SpilsBandCallback None
        | _ -> failwith "spils solver not in use"

      external get_work_space : serial_session -> int * int
        = "c_cvode_bandprec_get_work_space"

      external get_num_rhs_evals : serial_session -> int
        = "c_cvode_bandprec_get_num_rhs_evals"
    end
  end

module Alternate =
  struct
    include AlternateTypes

    external c_set_alternate
      : ('data, 'kind) session -> bool -> bool -> unit
      = "c_cvode_set_alternate"

    type gammas = {
      gamma : float;
      gammap : float;
    }

    external get_gammas : ('data, 'kind) session -> gammas
      = "c_cvode_get_gamma"

    let make_solver f s nv =
      let { linit; lsetup; lsolve } as cb = f s nv in
      c_set_alternate s (linit <> None) (lsetup <> None);
      s.ls_callbacks <- AlternateCallback cb

  end

external c_set_functional : ('a, 'k) session -> unit
  = "c_cvode_set_functional"

let set_iter_type session nv iter =
  session.ls_callbacks <- NoCallbacks;
  match iter with
  | Functional -> c_set_functional session
  | Newton linsolv -> linsolv session nv
    (* Iter type will be set to CV_NEWTON in the functions that set the linear
       solver.  *)

external sv_tolerances  : ('a, 'k) session -> float -> ('a, 'k) nvector -> unit
    = "c_cvode_sv_tolerances"
external ss_tolerances  : ('a, 'k) session -> float -> float -> unit
    = "c_cvode_ss_tolerances"
external wf_tolerances  : ('a, 'k) session -> unit
    = "c_cvode_wf_tolerances"

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

external c_session_finalize : ('a, 'kind) session -> unit
    = "c_cvode_session_finalize"

let session_finalize s =
  Dls.invalidate_callback s;
  c_session_finalize s

external c_init
    : ('a, 'k) session Weak.t -> lmm -> ('a, 'k) iter -> ('a, 'k) nvector
      -> float -> (cvode_mem * c_weak_ref * cvode_file)
    = "c_cvode_init"

let init lmm iter tol f ?(roots=no_roots) t0 y0 =
  let (nroots, roots) = roots in
  if nroots < 0 then
    raise (Invalid_argument "number of root functions is negative");
  let weakref = Weak.create 1 in
  let cvode_mem, backref, err_file = c_init weakref lmm iter y0 t0 in
  (* cvode_mem and backref have to be immediately captured in a session and
     associated with the finalizer before we do anything else.  *)
  let session = {
          cvode        = cvode_mem;
          backref      = backref;
          nroots       = nroots;
          err_file     = err_file;

          exn_temp     = None;

          rhsfn        = f;
          rootsfn      = roots;
          errh         = dummy_errh;
          errw         = dummy_errw;

          ls_callbacks = NoCallbacks;

          sensext      = NoSensExt;
        } in
  Gc.finalise session_finalize session;
  Weak.set weakref 0 (Some session);
  (* Now the session is safe to use.  If any of the following fails and raises
     an exception, the GC will take care of freeing cvode_mem and backref.  *)
  if nroots > 0 then
    c_root_init session nroots;
  set_iter_type session y0 iter;
  set_tolerances session tol;
  session

let nroots { nroots } = nroots

external c_reinit
    : ('a, 'k) session -> float -> ('a, 'k) nvector -> unit
    = "c_cvode_reinit"
let reinit session ?iter_type ?roots t0 y0 =
  Dls.invalidate_callback session;
  c_reinit session t0 y0;
  (match iter_type with
   | None -> ()
   | Some iter_type -> set_iter_type session y0 iter_type);
  (match roots with
   | None -> ()
   | Some roots -> root_init session roots)

external get_root_info  : ('a, 'k) session -> Sundials.Roots.t -> unit
    = "c_cvode_get_root_info"

external solve_normal : ('a, 'k) session -> float -> ('a, 'k) nvector
                              -> float * Sundials.solver_result
    = "c_cvode_solve_normal"

external solve_one_step : ('a, 'k) session -> float -> ('a, 'k) nvector
                              -> float * Sundials.solver_result
    = "c_cvode_solve_one_step"

external get_dky
    : ('a, 'k) session -> float -> int -> ('a, 'k) nvector -> unit
    = "c_cvode_get_dky"

external get_integrator_stats : ('a, 'k) session -> integrator_stats
    = "c_cvode_get_integrator_stats"

external get_work_space         : ('a, 'k) session -> int * int
    = "c_cvode_get_work_space"

external get_num_steps          : ('a, 'k) session -> int
    = "c_cvode_get_num_steps"

external get_num_rhs_evals      : ('a, 'k) session -> int
    = "c_cvode_get_num_rhs_evals"

external get_num_lin_solv_setups : ('a, 'k) session -> int
    = "c_cvode_get_num_lin_solv_setups"

external get_num_err_test_fails : ('a, 'k) session -> int
    = "c_cvode_get_num_err_test_fails"

external get_last_order         : ('a, 'k) session -> int
    = "c_cvode_get_last_order"

external get_current_order      : ('a, 'k) session -> int
    = "c_cvode_get_current_order"

external get_actual_init_step   : ('a, 'k) session -> float
    = "c_cvode_get_actual_init_step"

external get_last_step          : ('a, 'k) session -> float
    = "c_cvode_get_last_step"

external get_current_step       : ('a, 'k) session -> float
    = "c_cvode_get_current_step"

external get_current_time       : ('a, 'k) session -> float
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

external set_error_file : ('a, 'k) session -> string -> bool -> unit
    = "c_cvode_set_error_file"

external set_err_handler_fn  : ('a, 'k) session -> unit
    = "c_cvode_set_err_handler_fn"

let set_err_handler_fn s ferrh =
  s.errh <- ferrh;
  set_err_handler_fn s

external clear_err_handler_fn  : ('a, 'k) session -> unit
    = "c_cvode_clear_err_handler_fn"

let clear_err_handler_fn s =
  s.errh <- dummy_errh;
  clear_err_handler_fn s

external set_max_ord            : ('a, 'k) session -> int -> unit
    = "c_cvode_set_max_ord"
external set_max_num_steps      : ('a, 'k) session -> int -> unit
    = "c_cvode_set_max_num_steps"
external set_max_hnil_warns     : ('a, 'k) session -> int -> unit
    = "c_cvode_set_max_hnil_warns"
external set_stab_lim_det       : ('a, 'k) session -> bool -> unit
    = "c_cvode_set_stab_lim_det"
external set_init_step          : ('a, 'k) session -> float -> unit
    = "c_cvode_set_init_step"
external set_min_step           : ('a, 'k) session -> float -> unit
    = "c_cvode_set_min_step"
external set_max_step           : ('a, 'k) session -> float -> unit
    = "c_cvode_set_max_step"
external set_stop_time          : ('a, 'k) session -> float -> unit
    = "c_cvode_set_stop_time"
external set_max_err_test_fails : ('a, 'k) session -> int -> unit
    = "c_cvode_set_max_err_test_fails"
external set_max_nonlin_iters   : ('a, 'k) session -> int -> unit
    = "c_cvode_set_max_nonlin_iters"
external set_max_conv_fails     : ('a, 'k) session -> int -> unit
    = "c_cvode_set_max_conv_fails"
external set_nonlin_conv_coef   : ('a, 'k) session -> float -> unit
    = "c_cvode_set_nonlin_conv_coef"

external set_root_direction'   : ('a, 'k) session -> Sundials.RootDirs.t -> unit
    = "c_cvode_set_root_direction"

let set_root_direction s rda =
  set_root_direction' s (Sundials.RootDirs.copy (nroots s) rda)

let set_all_root_directions s rd =
  set_root_direction' s (Sundials.RootDirs.make (nroots s) rd)

external set_no_inactive_root_warn      : ('a, 'k) session -> unit
    = "c_cvode_set_no_inactive_root_warn"

external get_num_stab_lim_order_reds    : ('a, 'k) session -> int
    = "c_cvode_get_num_stab_lim_order_reds"

external get_tol_scale_factor           : ('a, 'k) session -> float
    = "c_cvode_get_tol_scale_factor"

external get_err_weights : ('a, 'k) session -> ('a, 'k) nvector -> unit
    = "c_cvode_get_err_weights"

external get_est_local_errors : ('a, 'k) session -> ('a, 'k) nvector -> unit
    = "c_cvode_get_est_local_errors"

external get_num_nonlin_solv_iters      : ('a, 'k) session -> int
    = "c_cvode_get_num_nonlin_solv_iters"

external get_num_nonlin_solv_conv_fails : ('a, 'k) session -> int
    = "c_cvode_get_num_nonlin_solv_conv_fails"

external get_nonlin_solv_stats          : ('a, 'k) session -> int * int
    = "c_cvode_get_nonlin_solv_stats"

external get_num_g_evals                : ('a, 'k) session -> int
    = "c_cvode_get_num_g_evals"


(* Callbacks *)

let call_rhsfn session t y y' =
  let session = read_weak_ref session in
  adjust_retcode session true (session.rhsfn t y) y'

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
                   "This exception will not be propagated: " ^
                   Printexc.to_string e)

let call_precsolvefn session jac r z =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | SpilsCallback { Spils.prec_solve_fn = f } ->
      adjust_retcode session true (f jac r) z
  | _ -> assert false

let call_precsetupfn session jac jok gamma =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | SpilsCallback { Spils.prec_setup_fn = Some f } ->
      adjust_retcode_and_bool session (f jac jok) gamma
  | _ -> assert false

let call_jactimesfn session jac v jv =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | SpilsCallback { Spils.jac_times_vec_fn = Some f }
  | SpilsBandCallback (Some f) ->
      adjust_retcode session true (f jac v) jv
  | _ -> assert false

let call_linit session =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | AlternateCallback { linit = Some f } ->
      adjust_retcode session false f session
  | _ -> assert false

let call_lsetup session convfail ypred fpred tmp =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | AlternateCallback { lsetup = Some f } ->
      adjust_retcode_and_bool session (f session convfail ypred fpred) tmp
  | _ -> assert false

let call_lsolve session b weight ycur fcur =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | AlternateCallback { lsolve = f } ->
      adjust_retcode session true (f session b weight ycur) fcur
  | _ -> assert false

(* Let C code know about some of the values in this module.  *)
type fcn = Fcn : 'a -> fcn
external c_init_module : fcn array -> exn array -> unit =
  "c_cvode_init_module"

let _ =
  c_init_module
    (* Functions must be listed in the same order as
       callback_index in cvode_ml.c.  *)
    [|Fcn call_rhsfn;
      Fcn call_errw;
      Fcn call_errh;

      Fcn call_precsolvefn;
      Fcn call_precsetupfn;
      Fcn call_jactimesfn;

      Fcn call_linit;
      Fcn call_lsetup;
      Fcn call_lsolve;
    |]

    (* Exceptions must be listed in the same order as
       cvode_exn_index.  *)
    [|IllInput;
      TooClose;
      TooMuchWork;
      TooMuchAccuracy;
      ErrFailure;
      ConvergenceFailure;
      LinearInitFailure;
      LinearSetupFailure;
      LinearSolveFailure;
      RhsFuncFailure;
      FirstRhsFuncFailure;
      RepeatedRhsFuncFailure;
      UnrecoverableRhsFuncFailure;
      RootFuncFailure;
      BadK;
      BadT;
      BadDky;
    |]
