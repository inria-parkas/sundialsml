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
include Ida_impl

(* "Simulate" Linear Solvers in Sundials < 3.0.0 *)
let in_compat_mode2 =
  match Config.sundials_version with
  | 2,_,_ -> true
  | _ -> false

(* "Simulate" Nonlinear Solvers in Sundials < 4.0.0 *)
let in_compat_mode2_3 =
  match Config.sundials_version with
  | 2,_,_ -> true
  | 3,_,_ -> true
  | _ -> false

(*
 * NB: The order of variant constructors and record fields is important!
 *     If these types are changed or augmented, the corresponding declarations
 *     in ida_ml.h (and code in ida_ml.c) must also be updated.
 *)

(* Solver exceptions *)
exception IllInput
exception TooMuchWork
exception TooMuchAccuracy
exception ErrFailure
exception ConvergenceFailure
exception LinearInitFailure
exception LinearSetupFailure of exn option
exception LinearSolveFailure of exn option
exception NonlinearInitFailure
exception NonlinearSetupFailure
exception NonlinearSetupRecoverable
exception ResFuncFailure
exception FirstResFuncFailure
exception RepeatedResFuncFailure
exception RootFuncFailure
exception ConstraintFailure

(* Initial condition calculator exceptions *)
exception LinesearchFailure
exception NoRecovery
exception BadEwt

(* get_dky exceptions *)
exception BadK
exception BadT

(* SetId exceptions *)
exception IdNotSet

exception VectorOpErr

let no_roots = (0, dummy_rootsfn)

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

module VarId = struct (* {{{ *)
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
end (* }}} *)

external c_root_init : ('a, 'k) session -> int -> unit
    = "sunml_ida_root_init"

let root_init session (nroots, rootsfn) =
  c_root_init session nroots;
  session.rootsfn <- rootsfn

module Dls = struct (* {{{ *)
  include DirectTypes
  include LinearSolver.Direct

  (* Sundials < 3.0.0 *)
  external c_dls_dense : 'k serial_session -> int -> bool -> unit
    = "sunml_ida_dls_dense"

  (* Sundials < 3.0.0 *)
  external c_dls_lapack_dense : 'k serial_session -> int -> bool -> unit
    = "sunml_ida_dls_lapack_dense"

  (* Sundials < 3.0.0 *)
  external c_dls_band : 'k serial_session -> int -> int -> int -> bool -> unit
    = "sunml_ida_dls_band"

  (* Sundials < 3.0.0 *)
  external c_dls_lapack_band : 'k serial_session -> int -> int -> int -> bool
                             -> unit
    = "sunml_ida_dls_lapack_band"

  (* Sundials < 3.0.0 *)
  external c_klu
    : 'k serial_session -> 's Matrix.Sparse.sformat -> int -> int -> unit
    = "sunml_ida_klu_init"

  (* Sundials < 3.0.0 *)
  external c_klu_set_ordering
    : 'k serial_session -> LinearSolver.Direct.Klu.ordering -> unit
    = "sunml_ida_klu_set_ordering"

  (* Sundials < 3.0.0 *)
  external c_klu_reinit : 'k serial_session -> int -> int -> unit
    = "sunml_ida_klu_reinit"

  (* Sundials < 3.0.0 *)
  external c_superlumt : 'k serial_session -> int -> int -> int -> unit
    = "sunml_ida_superlumt_init"

  (* Sundials < 3.0.0 *)
  external c_superlumt_set_ordering
    : 'k serial_session -> LinearSolver.Direct.Superlumt.ordering -> unit
    = "sunml_ida_superlumt_set_ordering"

  (* Sundials < 3.0.0 *)
  let klu_set_ordering session ordering =
    match session.ls_callbacks with
    | SlsKluCallback _ | BSlsKluCallback _ | BSlsKluCallbackSens _ ->
        c_klu_set_ordering session ordering
    | _ -> ()

  (* Sundials < 3.0.0 *)
  let klu_reinit session n onnz =
    match session.ls_callbacks with
    | SlsKluCallback _ | BSlsKluCallback _ | BSlsKluCallbackSens _ ->
        c_klu_reinit session n (match onnz with None -> 0 | Some nnz -> nnz)
    | _ -> ()

  (* Sundials < 3.0.0 *)
  let superlumt_set_ordering session ordering =
    match session.ls_callbacks with
    | SlsSuperlumtCallback _ | BSlsSuperlumtCallback _
    | BSlsSuperlumtCallbackSens _ ->
        c_superlumt_set_ordering session ordering
    | _ -> ()

  module LSD = LSI.Direct

  (* Sundials < 3.0.0 *)
  let make_compat (type s) (type tag) hasjac
        (solver : (s, 'nd, 'nk, tag) LSD.solver)
        (mat : ('k, s, 'nd, 'nk) Matrix.t) session =
    match solver with
    | LSD.Dense ->
        let m, n = Matrix.(Dense.size (unwrap mat)) in
        if m <> n then raise LinearSolver.MatrixNotSquare;
        c_dls_dense session m hasjac
    | LSD.LapackDense ->
        let m, n = Matrix.(Dense.size (unwrap mat)) in
        if m <> n then raise LinearSolver.MatrixNotSquare;
        c_dls_lapack_dense session m hasjac

    | LSD.Band ->
        let open Matrix.Band in
        let { n; mu; ml } = dims (Matrix.unwrap mat) in
        c_dls_band session n mu ml hasjac
    | LSD.LapackBand ->
        let open Matrix.Band in
        let { n; mu; ml } = dims (Matrix.unwrap mat) in
        c_dls_lapack_band session n mu ml hasjac

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
    = "sunml_ida_dls_set_linear_solver"

  (* 4.0.0 <= Sundials *)
  external c_set_linear_solver
    : ('d, 'k) session
      -> ('m, 'd, 'k) LSD.cptr
      -> ('mk, 'm, 'd, 'k) Matrix.t option
      -> bool
      -> unit
    = "sunml_ida_set_linear_solver"

  let solver ?jac ((LSD.S { LSD.rawptr; LSD.solver; LSD.matrix }) as ls)
             session nv =
    set_ls_callbacks ?jac solver matrix session;
    if in_compat_mode2
       then make_compat (jac <> None) solver matrix session
    else if in_compat_mode2_3
         then c_dls_set_linear_solver session rawptr matrix (jac <> None)
    else c_set_linear_solver session rawptr (Some matrix) (jac <> None);
    LSD.attach ls;
    session.ls_solver <- LSI.DirectSolver ls

  (* Sundials < 3.0.0 *)
  let invalidate_callback session =
    if in_compat_mode2 then
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
      = "sunml_ida_dls_get_work_space"

  let get_work_space s =
    if in_compat_mode2_3 then ls_check_direct s;
    get_work_space s

  external c_get_num_jac_evals : 'k serial_session -> int
      = "sunml_ida_get_num_jac_evals"

  external c_klu_get_num_jac_evals : 'k serial_session -> int
    = "sunml_ida_klu_get_num_jac_evals"

  external c_superlumt_get_num_jac_evals : 'k serial_session -> int
    = "sunml_ida_superlumt_get_num_jac_evals"

  let compat_get_num_jac_evals s =
    match s.ls_callbacks with
    | SlsKluCallback _ | BSlsKluCallback _ | BSlsKluCallbackSens _ ->
        c_klu_get_num_jac_evals s
    | SlsSuperlumtCallback _ | BSlsSuperlumtCallback _
    | BSlsSuperlumtCallbackSens _ -> c_superlumt_get_num_jac_evals s
    | _ -> c_get_num_jac_evals s

  let get_num_jac_evals s =
    if in_compat_mode2_3 then ls_check_direct s;
    if in_compat_mode2 then compat_get_num_jac_evals s else
    c_get_num_jac_evals s

  external get_num_lin_res_evals : 'k serial_session -> int
      = "sunml_ida_dls_get_num_lin_res_evals"

  let get_num_lin_res_evals s =
    if in_compat_mode2_3 then ls_check_direct s;
    get_num_lin_res_evals s
end (* }}} *)

module Spils = struct (* {{{ *)
  include SpilsTypes
  include LinearSolver.Iterative

  (* Sundials < 3.0.0 *)
  external c_spgmr
    : ('a, 'k) session -> int -> unit
    = "sunml_ida_spils_spgmr"

  (* Sundials < 3.0.0 *)
  external c_spbcgs
    : ('a, 'k) session -> int -> unit
    = "sunml_ida_spils_spbcgs"

  (* Sundials < 3.0.0 *)
  external c_sptfqmr
    : ('a, 'k) session -> int -> unit
    = "sunml_ida_spils_sptfqmr"

  (* Sundials < 3.0.0 *)
  external c_set_gs_type
    : ('a, 'k) session -> LSI.Iterative.gramschmidt_type -> unit
    = "sunml_ida_spils_set_gs_type"

  (* Sundials < 3.0.0 *)
  external c_set_maxl : ('a, 'k) session -> int -> unit
    = "sunml_ida_spils_set_maxl"

  (* Sundials < 3.0.0 *)
  external c_set_max_restarts : ('a, 'k) session -> int -> unit
    = "sunml_ida_spils_set_max_restarts"

  let old_set_maxl s maxl =
    ls_check_spils s;
    c_set_maxl s maxl

  let old_set_gs_type s t =
    ls_check_spils s;
    c_set_gs_type s t

  let old_set_max_restarts s t =
    ls_check_spils s;
    c_set_max_restarts s t

  external c_set_jac_times : ('a, 'k) session -> bool -> bool -> unit
    = "sunml_ida_set_jac_times"

  external c_set_preconditioner
    : ('a, 'k) session -> bool -> unit
    = "sunml_ida_set_preconditioner"

  external c_spils_set_linear_solver
    : ('a, 'k) session -> ('a, 'k) LSI.Iterative.cptr -> unit
    = "sunml_ida_spils_set_linear_solver"

  (* 4.0.0 <= Sundials *)
  external c_set_linear_solver
    : ('d, 'k) session
      -> ('d, 'k) LSI.Iterative.cptr
      -> ('mk, 'm, 'd, 'k) Matrix.t option
      -> bool
      -> unit
    = "sunml_ida_set_linear_solver"

  let init_preconditioner solve setup session nv =
    c_set_preconditioner session (setup <> None);
    session.ls_precfns <- PrecFns { prec_solve_fn = solve;
                                    prec_setup_fn = setup }

  let prec_none = LSI.Iterative.(PrecNone,
                    fun session nv -> session.ls_precfns <- NoPrecFns)

  let prec_left ?setup solve  = LSI.Iterative.(PrecLeft,
                                            init_preconditioner solve setup)

  let check_prec_type prec_type =
    let open LSI.Iterative in
    match prec_type with
    | PrecNone | PrecLeft -> true
    | PrecRight | PrecBoth -> false

  let solver (type s)
        ({ LSI.Iterative.rawptr;
           LSI.Iterative.solver;
           LSI.Iterative.compat; } as lsolver)
        ?jac_times_vec (prec_type, set_prec) session nv =
    let jac_times_setup, jac_times_vec =
      match jac_times_vec with None -> None, None
                             | Some (ojts, jtv) -> ojts, Some jtv in
    if in_compat_mode2 then begin
      if jac_times_setup <> None then
        raise Config.NotImplementedBySundialsVersion;
      let open LSI.Iterative in
      lsolver.check_prec_type <- check_prec_type;
      (match (solver : ('nd, 'nk, s) solver) with
       | Spgmr ->
           c_spgmr session compat.maxl;
           (match compat.gs_type with None -> ()
                                    | Some t -> c_set_gs_type session t);
           (match compat.max_restarts with None -> ()
                                      | Some t -> c_set_max_restarts session t);
           compat.set_gs_type <- old_set_gs_type session;
           compat.set_max_restarts <- old_set_max_restarts session
       | Spbcgs ->
           c_spbcgs session compat.maxl;
           compat.set_maxl <- old_set_maxl session
       | Sptfqmr ->
           c_sptfqmr session compat.maxl;
           compat.set_maxl <- old_set_maxl session
       | _ -> raise Config.NotImplementedBySundialsVersion);
      session.ls_solver <- LSI.IterativeSolver lsolver;
      set_prec session nv;
      session.ls_callbacks <- SpilsCallback (jac_times_vec, None);
      if jac_times_vec <> None then c_set_jac_times session true false
    end else
      if in_compat_mode2_3 then c_spils_set_linear_solver session rawptr
      else c_set_linear_solver session rawptr None false;
      LSI.Iterative.attach lsolver;
      session.ls_solver <- LSI.IterativeSolver lsolver;
      LSI.Iterative.(c_set_prec_type rawptr solver prec_type false);
      set_prec session nv;
      session.ls_callbacks <- SpilsCallback (jac_times_vec, jac_times_setup);
      if jac_times_setup <> None || jac_times_vec <> None then
        c_set_jac_times session (jac_times_setup <> None)
                                (jac_times_vec <> None)

  let set_jac_times s ?jac_times_setup f =
    if in_compat_mode2 && jac_times_setup <> None then
        raise Config.NotImplementedBySundialsVersion;
    match s.ls_callbacks with
    | SpilsCallback _ ->
        c_set_jac_times s (jac_times_setup <> None) true;
        s.ls_callbacks <- SpilsCallback (Some f, jac_times_setup)
    | _ -> raise LinearSolver.InvalidLinearSolver

  let clear_jac_times s =
    match s.ls_callbacks with
    | SpilsCallback _ ->
        c_set_jac_times s false false;
        s.ls_callbacks <- SpilsCallback (None, None)
    | _ -> raise LinearSolver.InvalidLinearSolver

  let set_preconditioner s ?setup solve =
    match s.ls_callbacks with
    | SpilsCallback _ ->
        c_set_preconditioner s (setup <> None);
        s.ls_precfns <- PrecFns { prec_setup_fn = setup;
                                  prec_solve_fn = solve }
    | _ -> raise LinearSolver.InvalidLinearSolver

  external set_eps_lin            : ('a, 'k) session -> float -> unit
      = "sunml_ida_set_eps_lin"

  let set_eps_lin s epsl =
    if in_compat_mode2_3 then ls_check_spils s;
    set_eps_lin s epsl

  external set_increment_factor   : ('a, 'k) session -> float -> unit
      = "sunml_ida_set_increment_factor"

  let set_increment_factor s dqincfac =
    if in_compat_mode2_3 then ls_check_spils s;
    set_increment_factor s dqincfac

  external get_num_lin_iters      : ('a, 'k) session -> int
      = "sunml_ida_get_num_lin_iters"

  let get_num_lin_iters s =
    if in_compat_mode2_3 then ls_check_spils s;
    get_num_lin_iters s

  external get_num_lin_conv_fails     : ('a, 'k) session -> int
      = "sunml_ida_get_num_lin_conv_fails"

  let get_num_lin_conv_fails s =
    if in_compat_mode2_3 then ls_check_spils s;
    get_num_lin_conv_fails s

  external get_work_space         : ('a, 'k) session -> int * int
      = "sunml_ida_spils_get_work_space"

  let get_work_space s =
    if in_compat_mode2_3 then ls_check_spils s;
    get_work_space s

  external get_num_prec_evals     : ('a, 'k) session -> int
      = "sunml_ida_get_num_prec_evals"

  let get_num_prec_evals s =
    if in_compat_mode2_3 then ls_check_spils s;
    get_num_prec_evals s

  external get_num_prec_solves    : ('a, 'k) session -> int
      = "sunml_ida_get_num_prec_solves"

  let get_num_prec_solves s =
    if in_compat_mode2_3 then ls_check_spils s;
    get_num_prec_solves s

  external get_num_jtsetup_evals   : ('a, 'k) session -> int
      = "sunml_ida_get_num_jtsetup_evals"

  let get_num_jtsetup_evals s =
    if in_compat_mode2_3 then ls_check_spils s;
    get_num_jtsetup_evals s

  external get_num_jtimes_evals   : ('a, 'k) session -> int
      = "sunml_ida_get_num_jtimes_evals"

  let get_num_jtimes_evals s =
    if in_compat_mode2_3 then ls_check_spils s;
    get_num_jtimes_evals s

  external get_num_lin_res_evals      : ('a, 'k) session -> int
      = "sunml_ida_spils_get_num_lin_res_evals"

  let get_num_lin_res_evals s =
    if in_compat_mode2_3 then ls_check_spils s;
    get_num_lin_res_evals s

end (* }}} *)

external sv_tolerances
    : ('a, 'k) session -> float -> ('a, 'k) Nvector.t -> unit
    = "sunml_ida_sv_tolerances"
external ss_tolerances  : ('a, 'k) session -> float -> float -> unit
    = "sunml_ida_ss_tolerances"
external wf_tolerances  : ('a, 'k) session -> unit
    = "sunml_ida_wf_tolerances"

type 'a error_fun = 'a -> 'a -> unit

type ('a, 'k) tolerance =
  | SStolerances of float * float
  | SVtolerances of float * ('a, 'k) Nvector.t
  | WFtolerances of 'a error_fun

let default_tolerances = SStolerances (1.0e-4, 1.0e-8)

let set_tolerances s tol =
  match tol with
  | SStolerances (rel, abs) -> (s.errw <- dummy_errw; ss_tolerances s rel abs)
  | SVtolerances (rel, abs) -> (if Sundials_configuration.safe then s.checkvec abs;
                                s.errw <- dummy_errw; sv_tolerances s rel abs)
  | WFtolerances ferrw -> (s.errw <- ferrw; wf_tolerances s)

external c_set_id : ('a,'k) session -> ('a,'k) Nvector.t -> unit
  = "sunml_ida_set_id"

let set_id s id =
  if Sundials_configuration.safe then s.checkvec id;
  c_set_id s id;
  s.id_set <- true

external c_session_finalize : ('a, 'kind) session -> unit
    = "sunml_ida_session_finalize"

let session_finalize s =
  Dls.invalidate_callback s;
  (match s.nls_solver with
   | None -> ()
   | Some nls -> NLSI.detach nls);
  (match s.sensext with
   | FwdSensExt { fnls_solver = Some nls } -> NLSI.detach nls
   | _ -> ());
  c_session_finalize s

(* Sundials >= 4.0.0 *)
external c_set_nonlinear_solver
    : ('d, 'k) session
      -> ('d, 'k, (('d, 'k) session) NLSI.integrator) NLSI.cptr
      -> unit
    = "sunml_ida_set_nonlinear_solver"

external c_init : ('a, 'k) session Weak.t -> float
                  -> ('a, 'k) Nvector.t -> ('a, 'k) Nvector.t
                  -> (ida_mem * c_weak_ref)
    = "sunml_ida_init"

let init tol ?nlsolver ~lsolver resfn ?varid ?(roots=no_roots) t0 y y' =
  let (nroots, rootsfn) = roots in
  let checkvec = Nvector.check y in
  if Sundials_configuration.safe then
    (checkvec y';
     if nroots < 0 then invalid_arg "number of root functions is negative");
  let weakref = Weak.create 1 in
  (if in_compat_mode2_3 then
    match nlsolver with
    | Some nls when NLSI.(get_type nls <> RootFind) -> raise IllInput
    | _ -> ());
  let ida_mem, backref = c_init weakref t0 y y' in
  (* ida_mem and backref have to be immediately captured in a session and
     associated with the finalizer before we do anything else.  *)
  let session = { ida        = ida_mem;
                  backref    = backref;
                  nroots     = nroots;
                  checkvec   = checkvec;

                  exn_temp   = None;

                  id_set     = false;
                  resfn      = resfn;
                  rootsfn    = rootsfn;
                  errh       = dummy_errh;
                  errw       = dummy_errw;

                  ls_solver  = LSI.NoSolver;
                  ls_callbacks = NoCallbacks;
                  ls_precfns = NoPrecFns;

                  nls_solver = None;

                  sensext    = NoSensExt;
                }
  in
  Gc.finalise session_finalize session;
  Weak.set weakref 0 (Some session);
  (* Now the session is safe to use.  If any of the following fails and raises
     an exception, the GC will take care of freeing ida_mem and backref.  *)
  (match varid with
     None -> ()
   | Some x -> set_id session x);
  if nroots > 0 then
    c_root_init session nroots;
  set_tolerances session tol;
  lsolver session y;
  (match nlsolver with
   | Some ({ NLSI.rawptr = nlcptr } as nls) when not in_compat_mode2_3 ->
       NLSI.attach nls;
       session.nls_solver <- Some nls;
       c_set_nonlinear_solver session nlcptr
   | _ -> ());
  session

let get_num_roots { nroots } = nroots

external c_reinit
    : ('a, 'k) session -> float -> ('a, 'k) Nvector.t
      -> ('a, 'k) Nvector.t -> unit
    = "sunml_ida_reinit"

let reinit session ?nlsolver ?lsolver ?roots t0 y0 y'0 =
  if Sundials_configuration.safe then
    (session.checkvec y0;
     session.checkvec y'0);
  Dls.invalidate_callback session;
  c_reinit session t0 y0 y'0;
  if in_compat_mode2_3 then
    match nlsolver with
    | Some nls when NLSI.(get_type nls <> RootFind) -> raise IllInput
    | _ -> ()
  else
    match nlsolver with
    | Some ({ NLSI.rawptr = nlcptr } as nls) ->
        (match session.nls_solver with
         | None -> () | Some old_nls -> NLSI.detach old_nls);
        NLSI.attach nls;
        session.nls_solver <- Some nls;
        c_set_nonlinear_solver session nlcptr
    | _ -> ();
  (match lsolver with
   | None -> ()
   | Some linsolv -> linsolv session y0);
  (match roots with
   | None -> ()
   | Some roots -> root_init session roots)

external get_root_info  : ('a, 'k) session -> Roots.t -> unit
    = "sunml_ida_get_root_info"

type solver_result =
  | Success             (** IDA_SUCCESS *)
  | RootsFound          (** IDA_ROOT_RETURN *)
  | StopTimeReached     (** IDA_TSTOP_RETURN *)

external c_solve_normal : ('a, 'k) session -> float
                          -> ('a, 'k) Nvector.t -> ('a,'k) Nvector.t
                          -> float * solver_result
    = "sunml_ida_solve_normal"

let solve_normal s t y yp =
  if Sundials_configuration.safe then
    (s.checkvec y;
     s.checkvec yp);
  c_solve_normal s t y yp

external c_solve_one_step : ('a, 'k) session -> float
                            -> ('a, 'k) Nvector.t-> ('a, 'k) Nvector.t
                            -> float * solver_result
    = "sunml_ida_solve_one_step"

let solve_one_step s t y yp =
  if Sundials_configuration.safe then
    (s.checkvec y;
     s.checkvec yp);
  c_solve_one_step s t y yp

external c_get_dky
    : ('a, 'k) session -> float -> int -> ('a, 'k) Nvector.t -> unit
    = "sunml_ida_get_dky"

let get_dky s y =
  if Sundials_configuration.safe then s.checkvec y;
  fun t k -> c_get_dky s t k y

external get_integrator_stats : ('a, 'k) session -> integrator_stats
    = "sunml_ida_get_integrator_stats"

external get_work_space         : ('a, 'k) session -> int * int
    = "sunml_ida_get_work_space"

external get_num_steps          : ('a, 'k) session -> int
    = "sunml_ida_get_num_steps"

external get_num_res_evals      : ('a, 'k) session -> int
    = "sunml_ida_get_num_res_evals"

external get_num_lin_solv_setups : ('a, 'k) session -> int
    = "sunml_ida_get_num_lin_solv_setups"

external get_num_err_test_fails : ('a, 'k) session -> int
    = "sunml_ida_get_num_err_test_fails"

external get_last_order         : ('a, 'k) session -> int
    = "sunml_ida_get_last_order"

external get_current_order      : ('a, 'k) session -> int
    = "sunml_ida_get_current_order"

external get_actual_init_step   : ('a, 'k) session -> float
    = "sunml_ida_get_actual_init_step"

external get_last_step          : ('a, 'k) session -> float
    = "sunml_ida_get_last_step"

external get_current_step       : ('a, 'k) session -> float
    = "sunml_ida_get_current_step"

external get_current_time       : ('a, 'k) session -> float
    = "sunml_ida_get_current_time"

let print_integrator_stats s oc =
  let stats = get_integrator_stats s
  in
    Printf.fprintf oc "num_steps = %d\n"           stats.num_steps;
    Printf.fprintf oc "num_res_evals = %d\n"       stats.num_res_evals;
    Printf.fprintf oc "num_lin_solv_setups = %d\n" stats.num_lin_solv_setups;
    Printf.fprintf oc "num_err_test_fails = %d\n"  stats.num_err_test_fails;
    Printf.fprintf oc "last_order = %d\n"          stats.last_order;
    Printf.fprintf oc "current_order = %d\n"       stats.current_order;
    Printf.fprintf oc "actual_init_step = %e\n"    stats.actual_init_step;
    Printf.fprintf oc "last_step = %e\n"           stats.last_step;
    Printf.fprintf oc "current_step = %e\n"        stats.current_step;
    Printf.fprintf oc "current_time = %e\n"        stats.current_time;

external set_error_file : ('a, 'k) session -> Logfile.t -> unit
    = "sunml_ida_set_error_file"

external set_err_handler_fn  : ('a, 'k) session -> unit
    = "sunml_ida_set_err_handler_fn"

let set_err_handler_fn s ferrh =
  s.errh <- ferrh;
  set_err_handler_fn s

external clear_err_handler_fn  : ('a, 'k) session -> unit
    = "sunml_ida_clear_err_handler_fn"

let clear_err_handler_fn s =
  s.errh <- dummy_errh;
  clear_err_handler_fn s

external set_max_ord            : ('a, 'k) session -> int -> unit
    = "sunml_ida_set_max_ord"
external set_max_num_steps      : ('a, 'k) session -> int -> unit
    = "sunml_ida_set_max_num_steps"
external set_init_step          : ('a, 'k) session -> float -> unit
    = "sunml_ida_set_init_step"
external set_max_step           : ('a, 'k) session -> float -> unit
    = "sunml_ida_set_max_step"
external set_stop_time          : ('a, 'k) session -> float -> unit
    = "sunml_ida_set_stop_time"
external set_max_err_test_fails : ('a, 'k) session -> int -> unit
    = "sunml_ida_set_max_err_test_fails"
external set_max_nonlin_iters   : ('a, 'k) session -> int -> unit
    = "sunml_ida_set_max_nonlin_iters"
external set_max_conv_fails     : ('a, 'k) session -> int -> unit
    = "sunml_ida_set_max_conv_fails"
external set_nonlin_conv_coef   : ('a, 'k) session -> float -> unit
    = "sunml_ida_set_nonlin_conv_coef"

external set_root_direction'   : ('a, 'k) session -> RootDirs.t -> unit
    = "sunml_ida_set_root_direction"

let set_root_direction s rda =
  set_root_direction' s (RootDirs.copy (get_num_roots s) rda)

let set_all_root_directions s rd =
  set_root_direction' s (RootDirs.make (get_num_roots s) rd)

external set_no_inactive_root_warn      : ('a, 'k) session -> unit
    = "sunml_ida_set_no_inactive_root_warn"
(*
   IDAGetNumStabLimOrderReds appears in the sundials 2.5.0 manual on
   p.52 but there's no such function in the implementation.  It's
   probably a leftover from earlier versions or something.

external get_num_stab_lim_order_reds    : ('a, 'k) session -> int
    = "c_ida_get_num_stab_lim_order_reds"
*)
external get_tol_scale_factor           : ('a, 'k) session -> float
    = "sunml_ida_get_tol_scale_factor"

external c_get_err_weights : ('a, 'k) session -> ('a, 'k) Nvector.t -> unit
    = "sunml_ida_get_err_weights"

let get_err_weights s ew =
  if Sundials_configuration.safe then s.checkvec ew;
  c_get_err_weights s ew

external c_get_est_local_errors : ('a, 'k) session -> ('a, 'k) Nvector.t -> unit
    = "sunml_ida_get_est_local_errors"

let get_est_local_errors s ew =
  if Sundials_configuration.safe then s.checkvec ew;
  c_get_est_local_errors s ew

external get_num_nonlin_solv_iters      : ('a, 'k) session -> int
    = "sunml_ida_get_num_nonlin_solv_iters"

external get_num_nonlin_solv_conv_fails : ('a, 'k) session -> int
    = "sunml_ida_get_num_nonlin_solv_conv_fails"

external get_nonlin_solv_stats          : ('a, 'k) session -> int * int
    = "sunml_ida_get_nonlin_solv_stats"

external get_num_g_evals                : ('a, 'k) session -> int
    = "sunml_ida_get_num_g_evals"

external c_set_constraints : ('a,'k) session -> ('a,'k) Nvector.t -> unit
  = "sunml_ida_set_constraints"

let set_constraints s nv =
  if Sundials_configuration.safe then s.checkvec nv;
  c_set_constraints s nv

external clear_constraints : ('a,'k) session -> unit
  = "sunml_ida_clear_constraints"

external c_set_suppress_alg : ('a,'k) session -> bool -> unit
  = "sunml_ida_set_suppress_alg"

let set_suppress_alg s ?varid v =
  (match varid with
   | None -> if v && not s.id_set then raise IdNotSet
   | Some x -> set_id s x);
  c_set_suppress_alg s v

external set_nonlin_conv_coef_ic : ('d, 'k) session -> float -> unit
  = "sunml_ida_set_nonlin_conv_coef_ic"

external set_max_num_steps_ic : ('d, 'k) session -> int -> unit
  = "sunml_ida_set_max_num_steps_ic"

external set_max_num_jacs_ic : ('d, 'k) session -> int -> unit
  = "sunml_ida_set_max_num_jacs_ic"

external set_max_num_iters_ic : ('d, 'k) session -> int -> unit
  = "sunml_ida_set_max_num_iters_ic"

external set_max_backs_ic : ('d, 'k) session -> int -> unit
  = "sunml_ida_set_max_backs_ic"

external set_line_search_ic : ('d, 'k) session -> bool -> unit
  = "sunml_ida_set_line_search_ic"

external set_step_tolerance_ic : ('d, 'k) session -> float -> unit
  = "sunml_ida_set_step_tolerance_ic"

external get_num_backtrack_ops : ('a,'k) session -> int
  = "sunml_ida_get_num_backtrack_ops"

external c_calc_ic_y : ('a,'k) session -> ('a,'k) Nvector.t option
                       -> float -> unit
  = "sunml_ida_calc_ic_y"

let calc_ic_y session ?y tout1 =
  if Sundials_configuration.safe then
    (match y with None -> () | Some x -> session.checkvec x);
  c_calc_ic_y session y tout1

external c_calc_ic_ya_yd' :
  ('a,'k) session -> ('a,'k) Nvector.t option -> ('a,'k) Nvector.t option
  -> float -> unit
  = "sunml_ida_calc_ic_ya_ydp"

let calc_ic_ya_yd' session ?y ?y' ?varid tout1 =
  if Sundials_configuration.safe then
    ((match y with None -> () | Some x -> session.checkvec x);
     (match y' with None -> () | Some x -> session.checkvec x));
  (match varid with
   | None -> if not session.id_set then raise IdNotSet
   | Some x -> set_id session x);
  c_calc_ic_ya_yd' session y y' tout1

(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "sunml_ida_init_module"

let _ =
  c_init_module
    (* Exceptions must be listed in the same order as
       ida_exn_index.  *)
    [|IllInput;
      TooMuchWork;
      TooMuchAccuracy;
      ErrFailure;
      ConvergenceFailure;
      LinearInitFailure;
      LinearSetupFailure None;
      LinearSolveFailure None;
      NonlinearInitFailure;
      NonlinearSetupFailure;
      NonlinearSetupRecoverable;
      ResFuncFailure;
      FirstResFuncFailure;
      RepeatedResFuncFailure;
      RootFuncFailure;
      ConstraintFailure;

      LinesearchFailure;
      NoRecovery;
      BadEwt;

      BadK;
      BadT;

      VectorOpErr;
    |]
