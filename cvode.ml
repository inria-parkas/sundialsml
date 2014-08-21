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
type serial_session = (real_array, Nvector_serial.kind) session

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
exception FirstRhsFuncErr
exception RepeatedRhsFuncErr
exception UnrecoverableRhsFuncErr
exception RootFuncFailure

(* get_dky exceptions *)
exception BadK
exception BadT
exception BadDky

let no_roots = (0, (fun _ _ _ -> ()))

type lmm =
  | Adams
  | BDF

type serial_linear_solver = (real_array, Nvector_serial.kind) linear_solver

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

exception StopTimeReached

let _ =
  List.iter (fun (nm, ex) -> Callback.register_exception nm ex)
  [
    ("cvode_StopTimeReached",         StopTimeReached);
    ("cvode_IllInput",                IllInput);
    ("cvode_TooClose",                TooClose);
    ("cvode_TooMuchWork",             TooMuchWork);
    ("cvode_TooMuchAccuracy",         TooMuchAccuracy);
    ("cvode_ErrFailure",              ErrFailure);
    ("cvode_ConvergenceFailure",      ConvergenceFailure);
    ("cvode_LinearInitFailure",       LinearInitFailure);
    ("cvode_LinearSetupFailure",      LinearSetupFailure);
    ("cvode_LinearSolveFailure",      LinearSolveFailure);
    ("cvode_RhsFuncFailure",          RhsFuncFailure);
    ("cvode_FirstRhsFuncErr",         FirstRhsFuncErr);
    ("cvode_RepeatedRhsFuncErr",      RepeatedRhsFuncErr);
    ("cvode_UnrecoverableRhsFuncErr", UnrecoverableRhsFuncErr);
    ("cvode_RootFuncFailure",         RootFuncFailure);

    ("cvode_BadK",                    BadK);
    ("cvode_BadT",                    BadT);
    ("cvode_BadDky",                  BadDky);
  ]

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

let call_precsolvefn session jac r z =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | SpilsCallback { prec_solve_fn = Some f } ->
      adjust_retcode session true (f jac r) z
  | _ -> assert false

let call_precsetupfn session jac jok gamma =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | SpilsCallback { prec_setup_fn = Some f } ->
      adjust_retcode_and_bool session (f jac jok) gamma
  | _ -> assert false

let call_jactimesfn session jac v jv =
  let session = read_weak_ref session in
  match session.ls_callbacks with
  | SpilsCallback { jac_times_vec_fn = Some f } ->
      adjust_retcode session true (f jac v) jv
  | _ -> assert false

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

let _ =
  Callback.register "c_cvode_call_rhsfn"         call_rhsfn;

  Callback.register "c_cvode_call_errh"          call_errh;
  Callback.register "c_cvode_call_errw"          call_errw;

  Callback.register "c_cvode_call_jacfn"         call_jacfn;
  Callback.register "c_cvode_call_bandjacfn"     call_bandjacfn;

  Callback.register "c_cvode_call_precsolvefn"   call_precsolvefn;
  Callback.register "c_cvode_call_precsetupfn"   call_precsetupfn;
  Callback.register "c_cvode_call_jactimesfn"    call_jactimesfn;

  Callback.register "c_cvode_call_linit"         call_linit;
  Callback.register "c_cvode_call_lsetup"        call_lsetup;
  Callback.register "c_cvode_call_lsolve"        call_lsolve;
  Callback.register "c_cvode_call_lfree"         call_lfree

external session_finalize : ('a, 'kind) session -> unit
    = "c_cvode_session_finalize"

external c_init
    : ('a, 'k) session Weak.t -> lmm -> ('a, 'k) iter -> ('a, 'k) nvector
      -> float -> (cvode_mem * c_weak_ref * cvode_file)
    = "c_cvode_init"

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
    type dense_jac_fn = (real_array triple_tmp, real_array) jacobian_arg
                                                  -> Dls.DenseMatrix.t -> unit

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
      let neqs = Sundials.RealArray.length (Sundials.unvec nv) in
      (session.ls_callbacks <- match jac with
                               | None -> NoCallbacks
                               | Some f -> DenseCallback f);
      c_dls_dense session neqs (jac <> None)

    let lapack_dense jac session nv =
      let neqs = Sundials.RealArray.length (Sundials.unvec nv) in
      (session.ls_callbacks <- match jac with
                               | None -> NoCallbacks
                               | Some f -> DenseCallback f);
      c_dls_lapack_dense session neqs (jac <> None)

    type band_jac_fn = bandrange
                        -> (real_array triple_tmp, real_array) jacobian_arg
                        -> Dls.BandMatrix.t -> unit

    let band p jac session nv =
      let neqs = Sundials.RealArray.length (Sundials.unvec nv) in
      (session.ls_callbacks <- match jac with
                              | None -> NoCallbacks
                              | Some f -> BandCallback f);
      c_dls_band (session, neqs) p.mupper p.mlower (jac <> None)

    let lapack_band p jac session nv =
      let neqs = Sundials.RealArray.length (Sundials.unvec nv) in
      (session.ls_callbacks <- match jac with
                               | None -> NoCallbacks
                               | Some f -> BandCallback f);
      c_dls_lapack_band (session, neqs) p.mupper p.mlower (jac <> None)

    let set_dense_jac_fn s fjacfn =
      s.ls_callbacks <- DenseCallback fjacfn;
      set_dense_jac_fn s

    external clear_dense_jac_fn : serial_session -> unit
        = "c_cvode_dls_clear_dense_jac_fn"

    let clear_dense_jac_fn s =
      match s.ls_callbacks with
      | DenseCallback _ -> (s.ls_callbacks <- NoCallbacks;
                            clear_dense_jac_fn s)
      | _ -> failwith "dense linear solver not in use"

    external set_band_jac_fn : serial_session -> unit
        = "c_cvode_dls_set_band_jac_fn"

    let set_band_jac_fn s fbandjacfn =
      s.ls_callbacks <- BandCallback fbandjacfn;
      set_band_jac_fn s

    external clear_band_jac_fn : serial_session -> unit
        = "c_cvode_dls_clear_band_jac_fn"

    let clear_band_jac_fn s =
      match s.ls_callbacks with
      | BandCallback _ -> (s.ls_callbacks <- NoCallbacks;
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
    type gramschmidt_type = Spils.gramschmidt_type =
      | ModifiedGS
      | ClassicalGS

    type preconditioning_type = Spils.preconditioning_type =
      | PrecNone
      | PrecLeft
      | PrecRight
      | PrecBoth

    type 'a solve_arg = 'a prec_solve_arg =
      {
        rhs : 'a;
        gamma : float;
        delta : float;
        left : bool;
      }

    type 'a callbacks = 'a spils_callbacks =
      {
        prec_solve_fn : (('a single_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                         -> 'a -> unit) option;

        prec_setup_fn : (('a triple_tmp, 'a) jacobian_arg -> bool
                            -> float -> bool) option;

        jac_times_vec_fn :
          (('a single_tmp, 'a) jacobian_arg
           -> 'a (* v *)
           -> 'a (* Jv *)
           -> unit) option;
      }

    let no_precond = { prec_solve_fn = None;
                       prec_setup_fn = None;
                       jac_times_vec_fn = None; }

    external c_spils_spgmr
      : ('a, 'k) session -> int -> Spils.preconditioning_type -> unit
      = "c_cvode_spils_spgmr"

    external c_spils_spbcg
      : ('a, 'k) session -> int -> Spils.preconditioning_type -> unit
      = "c_cvode_spils_spbcg"

    external c_spils_sptfqmr
      : ('a, 'k) session -> int -> Spils.preconditioning_type -> unit
      = "c_cvode_spils_sptfqmr"

    external c_spils_set_preconditioner
      : ('a, 'k) session -> bool -> bool -> unit
      = "c_cvode_spils_set_preconditioner"

    let set_precond session prec_type cb =
      match prec_type with
      | Spils.PrecNone -> ()
      | Spils.PrecLeft | Spils.PrecRight | Spils.PrecBoth ->
        match cb.prec_solve_fn with
        | None -> invalid_arg "preconditioning type is not PrecNone, but no \
                               solve function given"
        | Some solve_fn ->
            session.ls_callbacks <- SpilsCallback cb;
            c_spils_set_preconditioner session
              (cb.prec_setup_fn <> None)
              (cb.jac_times_vec_fn <> None)

    let spgmr maxl prec_type cb session nv =
      let maxl = match maxl with None -> 0 | Some ml -> ml in
      c_spils_spgmr session maxl prec_type;
      set_precond session prec_type cb
    
    let spbcg maxl prec_type cb session nv =
      let maxl = match maxl with None -> 0 | Some ml -> ml in
      c_spils_spbcg session maxl prec_type;
      set_precond session prec_type cb

    let sptfqmr maxl prec_type cb session nv =
      let maxl = match maxl with None -> 0 | Some ml -> ml in
      c_spils_sptfqmr session maxl prec_type;
      set_precond session prec_type cb

    let set_preconditioner s fprecsetupfn fprecsolvefn =
      (match s.ls_callbacks with
       | SpilsCallback cbs ->
           s.ls_callbacks <- SpilsCallback { cbs with
                               prec_setup_fn = fprecsetupfn;
                               prec_solve_fn = Some fprecsolvefn }
       | _ -> failwith "spils solver not in use");
      c_spils_set_preconditioner s (fprecsetupfn <> None) false

    external set_jac_times_vec_fn : ('a, 'k) session -> unit
        = "c_cvode_set_jac_times_vec_fn"

    let set_jac_times_vec_fn s fjactimesfn =
      (match s.ls_callbacks with
       | SpilsCallback cbs ->
           s.ls_callbacks <- SpilsCallback { cbs with
                               jac_times_vec_fn = Some fjactimesfn }
       | _ -> failwith "spils solver not in use");
      set_jac_times_vec_fn s

    external clear_jac_times_vec_fn : ('a, 'k) session -> unit
        = "c_cvode_clear_jac_times_vec_fn"

    let clear_jac_times_vec_fn s =
      (match s.ls_callbacks with
       | SpilsCallback cbs ->
           s.ls_callbacks <- SpilsCallback { cbs with jac_times_vec_fn = None }
       | _ -> failwith "spils solver not in use");
      clear_jac_times_vec_fn s

    external set_prec_type
        : ('a, 'k) session -> Spils.preconditioning_type -> unit
        = "c_cvode_spils_set_prec_type"

    external set_gs_type : ('a, 'k) session -> Spils.gramschmidt_type -> unit
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

    module Banded =
      struct
        external c_spils_banded_spgmr
          : (serial_session * int) -> int
                    -> int -> int -> Spils.preconditioning_type -> unit
          = "c_cvode_spils_banded_spgmr"

        external c_spils_banded_spbcg
          : (serial_session * int) -> int
                    -> int -> int -> Spils.preconditioning_type -> unit
          = "c_cvode_spils_banded_spbcg"

        external c_spils_banded_sptfqmr
          : (serial_session * int) -> int
                    -> int -> int -> Spils.preconditioning_type -> unit
          = "c_cvode_spils_banded_sptfqmr"

        let spgmr maxl prec_type br session nv =
          let neqs = Sundials.RealArray.length (Sundials.unvec nv) in
          let maxl = match maxl with None -> 0 | Some ml -> ml in
          c_spils_banded_spgmr (session, neqs) br.mupper br.mlower maxl prec_type

        let spbcg maxl prec_type br session nv =
          let neqs = Sundials.RealArray.length (Sundials.unvec nv) in
          let maxl = match maxl with None -> 0 | Some ml -> ml in
          c_spils_banded_spbcg (session, neqs) br.mupper br.mlower maxl prec_type

        let sptfqmr maxl prec_type br session nv =
          let neqs = Sundials.RealArray.length (Sundials.unvec nv) in
          let maxl = match maxl with None -> 0 | Some ml -> ml in
          c_spils_banded_sptfqmr (session, neqs) br.mupper br.mlower maxl prec_type

        external get_work_space : serial_session -> int * int
            = "c_cvode_bandprec_get_work_space"

        external get_num_rhs_evals : serial_session -> int
            = "c_cvode_bandprec_get_num_rhs_evals"
      end

  end

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
      = "c_cvode_set_alternate"

    let make_solver f s nv =
      let { linit; lsetup; lsolve; lfree } as cb = f s nv in
      c_set_alternate s (linit <> None) (lsetup <> None) (lfree <> None);
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

let init lmm iter tol f ?(roots=no_roots) ?(t0=0.) y0 =
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
  s.errh <- (fun _ -> ());
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
  set_root_direction' s (Sundials.RootDirs.copy_n (nroots s) rda)

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

