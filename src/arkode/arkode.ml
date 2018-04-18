(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2015 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)
include Arkode_impl

(* "Simulate" Linear Solvers in Sundials < 3.0.0 *)
let in_compat_mode =
  match Sundials.sundials_version with
  | 2,_,_ -> true
  | _ -> false

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
exception MassInitFailure
exception MassSetupFailure
exception MassSolveFailure
exception MassMultFailure
exception RhsFuncFailure
exception FirstRhsFuncFailure
exception RepeatedRhsFuncFailure
exception UnrecoverableRhsFuncFailure
exception RootFuncFailure
exception PostprocStepFailure

(* get_dky exceptions *)
exception BadK
exception BadT

let no_roots = (0, dummy_rootsfn)

type ('d, 'kind) iter =
  | Newton of ('d, 'kind) linear_solver
  | FixedPoint of int

type linearity =
  | Linear of bool
  | Nonlinear

type ('d, 'k) imex = {
    implicit: 'd rhsfn * ('d, 'k) iter * linearity;
    explicit: 'd rhsfn
  }

type ('d, 'k) problem =
  | Implicit of 'd rhsfn * ('d, 'k) iter * linearity
  | Explicit of 'd rhsfn
  | ImEx of ('d, 'k) imex

type integrator_stats = {
    num_steps           : int;
    exp_steps           : int;
    acc_steps           : int;
    step_attempts       : int;
    num_nfe_evals       : int;
    num_nfi_evals       : int;
    num_lin_solv_setups : int;
    num_err_test_fails  : int;
    actual_init_step    : float;
    last_step           : float;
    current_step        : float;
    current_time        : float
  }

external c_root_init : ('a, 'k) session -> int -> unit
    = "c_arkode_root_init"

let root_init session (nroots, rootsfn) =
  c_root_init session nroots;
  session.rootsfn <- rootsfn

module Dls = struct (* {{{ *)
  include DirectTypes

  (* Sundials < 3.0.0 *)
  external c_dls_dense : 'k serial_session -> int -> bool -> unit
    = "c_arkode_dls_dense"

  (* Sundials < 3.0.0 *)
  external c_dls_lapack_dense : 'k serial_session -> int -> bool -> unit
    = "c_arkode_dls_lapack_dense"

  (* Sundials < 3.0.0 *)
  external c_dls_band : 'k serial_session -> int -> int -> int -> bool -> unit
    = "c_arkode_dls_band"

  (* Sundials < 3.0.0 *)
  external c_dls_lapack_band
    : 'k serial_session -> int -> int -> int -> bool -> unit
    = "c_arkode_dls_lapack_band"

  (* Sundials < 3.0.0 *)
  external c_klu
    : 'k serial_session -> 's Matrix.Sparse.sformat -> int -> int -> unit
    = "c_arkode_klu_init"

  (* Sundials < 3.0.0 *)
  external c_klu_set_ordering
    : 'k serial_session -> Lsolver.Direct.Klu.ordering -> unit
    = "c_arkode_klu_set_ordering"

  (* Sundials < 3.0.0 *)
  external c_klu_reinit : 'k serial_session -> int -> int -> unit
    = "c_arkode_klu_reinit"

  (* Sundials < 3.0.0 *)
  external c_superlumt : 'k serial_session -> int -> int -> int -> unit
    = "c_arkode_superlumt_init"

  (* Sundials < 3.0.0 *)
  external c_superlumt_set_ordering
    : 'k serial_session -> Lsolver.Direct.Superlumt.ordering -> unit
    = "c_arkode_superlumt_set_ordering"

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
    | SlsSuperlumtCallback _ -> c_superlumt_set_ordering session ordering
    | _ -> ()

  (* Sundials < 3.0.0 *)
  let make_compat (type s) (type tag) hasjac
        (solver : (s, 'nd, 'nk, tag) Lsolver_impl.Direct.solver)
        (mat : ('k, s, 'nd, 'nk) Matrix.t) session =
    match solver with
    | Lsolver_impl.Direct.Dense ->
        let m, n = Matrix.(Dense.size (unwrap mat)) in
        if m <> n then raise Lsolver.MatrixNotSquare;
        c_dls_dense session m hasjac
    | Lsolver_impl.Direct.LapackDense ->
        let m, n = Matrix.(Dense.size (unwrap mat)) in
        if m <> n then raise Lsolver.MatrixNotSquare;
        c_dls_lapack_dense session m hasjac

    | Lsolver_impl.Direct.Band ->
        let open Matrix.Band in
        let { n; mu; ml } = dims (Matrix.unwrap mat) in
        c_dls_band session n mu ml hasjac
    | Lsolver_impl.Direct.LapackBand ->
        let open Matrix.Band in
        let { n; mu; ml } = dims (Matrix.unwrap mat) in
        c_dls_lapack_band session n mu ml hasjac

    | Lsolver_impl.Direct.Klu sinfo ->
        if not Sundials_config.klu_enabled
          then raise Sundials.NotImplementedBySundialsVersion;
        let smat = Matrix.unwrap mat in
        let m, n = Matrix.Sparse.size smat in
        let nnz, _ = Matrix.Sparse.dims smat in
        if m <> n then raise Lsolver.MatrixNotSquare;
        let open Lsolver_impl.Klu in
        sinfo.set_ordering <- klu_set_ordering session;
        sinfo.reinit <- klu_reinit session;
        c_klu session (Matrix.Sparse.sformat smat) m nnz;
        (match sinfo.ordering with None -> ()
                                 | Some o -> c_klu_set_ordering session o)

    | Lsolver_impl.Direct.Superlumt sinfo ->
        if not Sundials_config.superlumt_enabled
          then raise Sundials.NotImplementedBySundialsVersion;
        let smat = Matrix.unwrap mat in
        let m, n = Matrix.Sparse.size smat in
        let nnz, _ = Matrix.Sparse.dims smat in
        if m <> n then raise Lsolver.MatrixNotSquare;
        let open Lsolver_impl.Superlumt in
        sinfo.set_ordering <- superlumt_set_ordering session;
        c_superlumt session m nnz sinfo.num_threads;
        (match sinfo.ordering with None -> ()
                                 | Some o -> c_superlumt_set_ordering session o)

    | Lsolver_impl.Direct.Custom _ ->
        assert false

  let check_dqjac jac mat = let open Matrix in
    match get_id mat with
    | Dense | Band -> ()
    | _ -> if jac = None then invalid_arg "A Jacobian function is required"

  let set_ls_callbacks (type m) (type tag)
        ?jac (solver : (m, 'nd, 'nk, tag) Lsolver_impl.Direct.solver)
        (mat : ('mk, m, 'nd, 'nk) Matrix.t) session =
    let cb = { jacfn = (match jac with None -> no_callback | Some f -> f);
               jmat  = (None : m option) } in
    let open Lsolver_impl.Direct in
    begin match solver with
    | Dense ->
        session.ls_callbacks <- DlsDenseCallback (cb, Matrix.unwrap mat)
    | LapackDense ->
        session.ls_callbacks <- DlsDenseCallback (cb, Matrix.unwrap mat)
    | Band ->
        session.ls_callbacks <- DlsBandCallback (cb, Matrix.unwrap mat)
    | LapackBand ->
        session.ls_callbacks <- DlsBandCallback (cb, Matrix.unwrap mat)
    | Klu _ ->
        if jac = None then invalid_arg "Klu requires Jacobian function";
        session.ls_callbacks <- SlsKluCallback (cb, Matrix.unwrap mat)
    | Superlumt _ ->
        if jac = None then invalid_arg "Superlumt requires Jacobian function";
        session.ls_callbacks <- SlsSuperlumtCallback (cb, Matrix.unwrap mat)
    | Custom _ ->
        check_dqjac jac mat;
        session.ls_callbacks <- DirectCustomCallback (cb, Matrix.unwrap mat)
    end;
    session.ls_precfns <- NoPrecFns

  (* Sundials >= 3.0.0 *)
  external c_dls_set_linear_solver
    : 'k serial_session
      -> ('m, Nvector_serial.data, 'k) Lsolver_impl.Direct.cptr
      -> ('mk, 'm, Nvector_serial.data, 'k) Matrix.t
      -> bool
      -> unit
    = "c_arkode_dls_set_linear_solver"

  let make ({ Lsolver_impl.Direct.rawptr; Lsolver_impl.Direct.solver } as ls)
           ?jac mat session nv =
    set_ls_callbacks ?jac solver mat session;
    if in_compat_mode then make_compat (jac <> None) solver mat session
    else c_dls_set_linear_solver session rawptr mat (jac <> None);
    Lsolver_impl.Direct.attach ls;
    session.ls_solver <- Lsolver_impl.DirectSolver ls

  (* Sundials < 3.0.0 *)
  let invalidate_callback session =
    if in_compat_mode then
      match session.ls_callbacks with
      | DlsDenseCallback ({ jmat = Some d } as cb, _) ->
          Matrix.Dense.invalidate d;
          cb.jmat <- None
      | DlsBandCallback  ({ jmat = Some d } as cb, _) ->
          Matrix.Band.invalidate d;
          cb.jmat <- None
      | SlsKluCallback ({ jmat = Some d } as cb, _) ->
          Matrix.Sparse.invalidate d;
          cb.jmat <- None
      | SlsSuperlumtCallback ({ jmat = Some d } as cb, _) ->
          Matrix.Sparse.invalidate d;
          cb.jmat <- None
      | _ -> ()

  external get_work_space : 'k serial_session -> int * int
      = "c_arkode_dls_get_work_space"

  let get_work_space s =
    ls_check_direct s;
    get_work_space s

  external c_get_num_jac_evals : 'k serial_session -> int
    = "c_arkode_dls_get_num_jac_evals"

  (* Sundials < 3.0.0 *)
  external c_klu_get_num_jac_evals : 'k serial_session -> int
    = "c_arkode_klu_get_num_jac_evals"

  (* Sundials < 3.0.0 *)
  external c_superlumt_get_num_jac_evals : 'k serial_session -> int
    = "c_arkode_superlumt_get_num_jac_evals"

  let compat_get_num_jac_evals s =
    match s.ls_callbacks with
    | SlsKluCallback _ -> c_klu_get_num_jac_evals s
    | SlsSuperlumtCallback _ -> c_superlumt_get_num_jac_evals s
    | _ -> c_get_num_jac_evals s

  let get_num_jac_evals s =
    ls_check_direct s;
    if in_compat_mode then compat_get_num_jac_evals s else
    c_get_num_jac_evals s

  external get_num_rhs_evals : 'k serial_session -> int
      = "c_arkode_dls_get_num_rhs_evals"

  let get_num_rhs_evals s =
    ls_check_direct s;
    get_num_rhs_evals s

end (* }}} *)

module Spils = struct (* {{{ *)
  include SpilsTypes

  (* Sundials < 3.0.0 *)
  external c_spgmr
    : ('a, 'k) session
      -> int -> Lsolver_impl.Iterative.preconditioning_type -> unit
    = "c_arkode_spils_spgmr"

  (* Sundials < 3.0.0 *)
  external c_spbcgs
    : ('a, 'k) session
      -> int -> Lsolver_impl.Iterative.preconditioning_type -> unit
    = "c_arkode_spils_spbcgs"

  (* Sundials < 3.0.0 *)
  external c_sptfqmr
    : ('a, 'k) session
      -> int -> Lsolver_impl.Iterative.preconditioning_type -> unit
    = "c_arkode_spils_sptfqmr"

  (* Sundials < 3.0.0 *)
  external c_spfgmr
    : ('a, 'k) session
      -> int -> Lsolver_impl.Iterative.preconditioning_type -> unit
    = "c_arkode_spils_spfgmr"

  (* Sundials < 3.0.0 *)
  external c_pcg
    : ('a, 'k) session
      -> int -> Lsolver_impl.Iterative.preconditioning_type -> unit
    = "c_arkode_spils_pcg"

  (* Sundials < 3.0.0 *)
  external c_set_gs_type
    : ('a, 'k) session -> Lsolver_impl.Iterative.gramschmidt_type -> unit
    = "c_arkode_spils_set_gs_type"

  (* Sundials < 3.0.0 *)
  external set_maxl
    : ('a, 'k) session -> int -> unit
    = "c_arkode_spils_set_maxl"

  (* Sundials < 3.0.0 *)
  external c_set_prec_type
    : ('a, 'k) session -> Lsolver_impl.Iterative.preconditioning_type -> unit
    = "c_arkode_spils_set_prec_type"

  let old_set_maxl s maxl =
    ls_check_spils s;
    set_maxl s maxl

  let old_set_prec_type s t =
    ls_check_spils s;
    c_set_prec_type s t

  let old_set_gs_type s t =
    ls_check_spils s;
    c_set_gs_type s t

  external c_set_jac_times : ('a, 'k) session -> bool -> bool -> unit
    = "c_arkode_spils_set_jac_times"

  external c_set_preconditioner
    : ('a, 'k) session -> bool -> unit
    = "c_arkode_spils_set_preconditioner"

  external c_spils_set_linear_solver
    : ('a, 'k) session -> ('a, 'k) Lsolver_impl.Iterative.cptr -> unit
    = "c_cvode_spils_set_linear_solver"

  let init_preconditioner solve setup session nv =
    c_set_preconditioner session (setup <> None);
    session.ls_precfns <- PrecFns { prec_solve_fn = solve;
                                    prec_setup_fn = setup }

  let prec_none = Lsolver_impl.Iterative.(PrecNone,
                    fun session nv -> session.ls_precfns <- NoPrecFns)

  let prec_left ?setup solve  = Lsolver_impl.Iterative.(PrecLeft,
                                            init_preconditioner solve setup)

  let prec_right ?setup solve = Lsolver_impl.Iterative.(PrecRight,
                                            init_preconditioner solve setup)

  let prec_both ?setup solve  = Lsolver_impl.Iterative.(PrecBoth,
                                            init_preconditioner solve setup)

  let make (type s)
        ({ Lsolver_impl.Iterative.rawptr;
           Lsolver_impl.Iterative.solver;
           Lsolver_impl.Iterative.compat =
             ({ Lsolver_impl.Iterative.maxl;
                Lsolver_impl.Iterative.gs_type } as compat) } as ls)
        ?jac_times_vec (prec_type, set_prec) session nv =
    let jac_times_setup, jac_times_vec =
      match jac_times_vec with None -> None, None
                             | Some (ojts, jtv) -> ojts, Some jtv in
    if in_compat_mode then begin
      if jac_times_setup <> None then
        raise Sundials.NotImplementedBySundialsVersion;
      let open Lsolver_impl.Iterative in
      (match (solver : ('nd, 'nk, s) solver) with
       | Spgmr ->
           c_spgmr session maxl prec_type;
           (match gs_type with None -> () | Some t -> c_set_gs_type session t);
           compat.set_gs_type <- old_set_gs_type session;
           compat.set_prec_type <- old_set_prec_type session
       | Spfgmr ->
           c_spfgmr session maxl prec_type;
           (match gs_type with None -> () | Some t -> c_set_gs_type session t);
           compat.set_gs_type <- old_set_gs_type session;
           compat.set_prec_type <- old_set_prec_type session
       | Spbcgs ->
           c_spbcgs session maxl prec_type;
           compat.set_maxl <- old_set_maxl session;
           compat.set_prec_type <- old_set_prec_type session
       | Sptfqmr ->
           c_sptfqmr session maxl prec_type;
           compat.set_maxl <- old_set_maxl session;
           compat.set_prec_type <- old_set_prec_type session
       | Pcg ->
           c_pcg session maxl prec_type;
           compat.set_maxl <- old_set_maxl session;
           compat.set_prec_type <- old_set_prec_type session
       | Custom _ -> assert false);
      session.ls_solver <- Lsolver_impl.IterativeSolver ls;
      set_prec session nv;
      session.ls_callbacks <- SpilsCallback (jac_times_vec, None);
      if jac_times_vec <> None then c_set_jac_times session true false
    end else
      c_spils_set_linear_solver session rawptr;
      Lsolver_impl.Iterative.attach ls;
      session.ls_solver <- Lsolver_impl.IterativeSolver ls;
      Lsolver_impl.Iterative.(c_set_prec_type rawptr solver prec_type);
      set_prec session nv;
      session.ls_callbacks <- SpilsCallback (jac_times_vec, jac_times_setup);
      if jac_times_setup <> None || jac_times_vec <> None then
        c_set_jac_times session (jac_times_setup <> None)
                                (jac_times_vec <> None)

  let set_jac_times s ?jac_times_setup f =
    if in_compat_mode && jac_times_setup <> None then
        raise Sundials.NotImplementedBySundialsVersion;
    match s.ls_callbacks with
    | SpilsCallback _ ->
        c_set_jac_times s (jac_times_setup <> None) true;
        s.ls_callbacks <- SpilsCallback (Some f, jac_times_setup)
    | _ -> raise Sundials.InvalidLinearSolver

  let clear_jac_times s =
    match s.ls_callbacks with
    | SpilsCallback _ ->
        c_set_jac_times s false false;
        s.ls_callbacks <- SpilsCallback (None, None)
    | _ -> raise Sundials.InvalidLinearSolver

  let set_preconditioner s ?setup solve =
    match s.ls_callbacks with
    | SpilsCallback _ ->
        c_set_preconditioner s (setup <> None);
        s.ls_precfns <- PrecFns { prec_setup_fn = setup;
                                  prec_solve_fn = solve }
    | _ -> raise Sundials.InvalidLinearSolver

  external set_eps_lin            : ('a, 'k) session -> float -> unit
    = "c_arkode_spils_set_eps_lin"

  let set_eps_lin s epsl =
    ls_check_spils s;
    set_eps_lin s epsl

  external get_num_lin_iters      : ('a, 'k) session -> int
    = "c_arkode_spils_get_num_lin_iters"

  let get_num_lin_iters s =
    ls_check_spils s;
    get_num_lin_iters s

  external get_num_conv_fails     : ('a, 'k) session -> int
    = "c_arkode_spils_get_num_conv_fails"

  let get_num_conv_fails s =
    ls_check_spils s;
    get_num_conv_fails s

  external get_work_space         : ('a, 'k) session -> int * int
    = "c_arkode_spils_get_work_space"

  let get_work_space s =
    ls_check_spils s;
    get_work_space s

  external get_num_prec_evals     : ('a, 'k) session -> int
    = "c_arkode_spils_get_num_prec_evals"

  let get_num_prec_evals s =
    ls_check_spils s;
    get_num_prec_evals s

  external get_num_prec_solves    : ('a, 'k) session -> int
    = "c_arkode_spils_get_num_prec_solves"

  let get_num_prec_solves s =
    ls_check_spils s;
    get_num_prec_solves s

  external get_num_jtsetup_evals   : ('a, 'k) session -> int
    = "c_arkode_spils_get_num_jtsetup_evals"

  let get_num_jtsetup_evals s =
    ls_check_spils s;
    get_num_jtsetup_evals s

  external get_num_jtimes_evals   : ('a, 'k) session -> int
    = "c_arkode_spils_get_num_jtimes_evals"

  let get_num_jtimes_evals s =
    ls_check_spils s;
    get_num_jtimes_evals s

  external get_num_rhs_evals      : ('a, 'k) session -> int
    = "c_arkode_spils_get_num_rhs_evals"

  let get_num_rhs_evals s =
    ls_check_spils s;
    get_num_rhs_evals s

  module Banded = struct (* {{{ *)
    (* These fields are accessed from arkode_ml.c *)
    type bandrange = { mupper : int; mlower : int; }

    external c_set_preconditioner
      : ('a, 'k) session -> int -> int -> int -> unit
      = "c_arkode_spils_set_banded_preconditioner"
    (* Note: ARKBandPrecInit seems to be designed only to be called on
       a fresh spils solver (i.e. right after ARKSpgmr, ARKSpbcg,
       ARKSptfqmr, ARKSpfgmr, or ARKPcg).

       As of Sundials 2.5.0, calling

         ARKSpgmr -> ARKBandPrecInit -> ARKBandPrecInit

       triggers a memory leak.  Calling ARKodeReInit does NOT help.
       The only way to prevent leakage is to allocate a fresh spils
       instance, thus:

         ARKSpgmr -> ARKBandPrecInit -> ARKSpgmr -> ARKBandPrecInit.

       If you call

         ARKSpgmr -> ARKSpilsSetPreconditioner -> ARKBandPrecInit,

       nothing grave happens, but the memory associated with
       ARKBandPrecInit won't be freed until the spils solver is torn
       down.  If you call BandPrecInit -> SetPreconditioner ->
       BandPrecInit, you also get a memory leak.

       (Perhaps set_preconditioner should be hidden too?  In that
        case, we should somehow strip set_prec_type of the ability
        to change PrecNone to anything else.)

       set_jac_times_vec_fn, as well as similar functions for Dls
       solvers, are kept because they accept NULL to remove the
       callback.  This design clearly anticipates being called
       multiple times on the same solver instance.  *)

    let init_preconditioner bandrange session nv =
      let n = RealArray.length (Nvector.unwrap nv) in
      c_set_preconditioner session n bandrange.mupper bandrange.mlower;
      session.ls_precfns <- BandedPrecFns

    let prec_none =
      Lsolver_impl.Iterative.(PrecNone, fun session nv ->
                                          session.ls_precfns <- BandedPrecFns)
    let prec_left bandrange =
      Lsolver_impl.Iterative.(PrecLeft,  init_preconditioner bandrange)
    let prec_right bandrange =
      Lsolver_impl.Iterative.(PrecRight, init_preconditioner bandrange)
    let prec_both bandrange =
      Lsolver_impl.Iterative.(PrecBoth,  init_preconditioner bandrange)

    external get_work_space : 'k serial_session -> int * int
      = "c_arkode_bandprec_get_work_space"

    let get_work_space s =
      ls_check_spils_band s;
      get_work_space s

    external get_num_rhs_evals : 'k serial_session -> int
      = "c_arkode_bandprec_get_num_rhs_evals"

    let get_num_rhs_evals s =
      ls_check_spils_band s;
      get_num_rhs_evals s
  end (* }}} *)

end (* }}} *)

module Alternate = struct (* {{{ *)
  include AlternateTypes

  external c_set_alternate
    : ('data, 'kind) session -> bool -> bool -> unit
    = "c_arkode_set_alternate"

  type gammas = {
    gamma : float;
    gammap : float;
  }

  external get_gammas : ('data, 'kind) session -> gammas
    = "c_arkode_get_gamma"

  let make f s nv =
    let { linit; lsetup; lsolve } as cb = f s nv in
    c_set_alternate s (linit <> None) (lsetup <> None);
    s.ls_precfns <- NoPrecFns;
    s.ls_callbacks <- AlternateCallback cb

end (* }}} *)

module Mass = struct (* {{{ *)
  include MassTypes

  module Dls = struct (* {{{ *)
    include MassTypes.Direct'

    (* Sundials < 3.0.0 *)
    external c_dls_mass_dense : 'k serial_session -> int -> unit
      = "c_arkode_dls_mass_dense"

    (* Sundials < 3.0.0 *)
    external c_dls_mass_lapack_dense : 'k serial_session -> int -> unit
      = "c_arkode_dls_mass_lapack_dense"

    (* Sundials < 3.0.0 *)
    external c_dls_mass_band : 'k serial_session -> int -> int -> int -> unit
      = "c_arkode_dls_mass_band"

    (* Sundials < 3.0.0 *)
    external c_dls_mass_lapack_band
      : 'k serial_session -> int -> int -> int -> unit
      = "c_arkode_dls_mass_lapack_band"

    (* Sundials < 3.0.0 *)
    external c_mass_klu
      : 'k serial_session -> 's Matrix.Sparse.sformat -> int -> int -> unit
      = "c_arkode_mass_klu_init"

    (* Sundials < 3.0.0 *)
    external c_klu_set_ordering
      : 'k serial_session -> Lsolver_impl.Klu.ordering -> unit
      = "c_arkode_mass_klu_set_ordering"

    (* Sundials < 3.0.0 *)
    external c_klu_reinit : 'k serial_session -> int -> int -> unit
      = "c_arkode_mass_klu_reinit"

    (* Sundials < 3.0.0 *)
    external c_mass_superlumt : 'k serial_session -> int -> int -> int -> unit
      = "c_arkode_mass_superlumt_init"

    (* Sundials < 3.0.0 *)
    external c_superlumt_set_ordering
      : 'k serial_session -> Lsolver_impl.Superlumt.ordering -> unit
      = "c_arkode_mass_superlumt_set_ordering"

    (* Sundials < 3.0.0 *)
    let klu_set_ordering session ordering =
      match session.mass_callbacks with
      | SlsKluMassCallback _ -> c_klu_set_ordering session ordering
      | _ -> ()

    (* Sundials < 3.0.0 *)
    let klu_reinit session n onnz =
      match session.mass_callbacks with
      | SlsKluMassCallback _ ->
          c_klu_reinit session n (match onnz with None -> 0 | Some nnz -> nnz)
      | _ -> ()

    (* Sundials < 3.0.0 *)
    let superlumt_set_ordering session ordering =
      match session.mass_callbacks with
      | SlsSuperlumtMassCallback _ -> c_superlumt_set_ordering session ordering
      | _ -> ()

    (* Sundials < 3.0.0 *)
    let make_compat (type s) (type tag)
          (solver : (s, 'nd, 'nk, tag) Lsolver_impl.Direct.solver)
          (mat : ('k, s, 'nd, 'nk) Matrix.t) session =
      match solver with
      | Lsolver_impl.Direct.Dense ->
          let m, n = Matrix.(Dense.size (unwrap mat)) in
          if m <> n then raise Lsolver.MatrixNotSquare;
          c_dls_mass_dense session m
      | Lsolver_impl.Direct.LapackDense ->
          let m, n = Matrix.(Dense.size (unwrap mat)) in
          if m <> n then raise Lsolver.MatrixNotSquare;
          c_dls_mass_lapack_dense session m

      | Lsolver_impl.Direct.Band ->
          let open Matrix.Band in
          let { n; mu; ml } = dims (Matrix.unwrap mat) in
          c_dls_mass_band session n mu ml
      | Lsolver_impl.Direct.LapackBand ->
          let open Matrix.Band in
          let { n; mu; ml } = dims (Matrix.unwrap mat) in
          c_dls_mass_lapack_band session n mu ml

      | Lsolver_impl.Direct.Klu sinfo ->
          if not Sundials_config.klu_enabled
            then raise Sundials.NotImplementedBySundialsVersion;
          let smat = Matrix.unwrap mat in
          let m, n = Matrix.Sparse.size smat in
          let nnz, _ = Matrix.Sparse.dims smat in
          if m <> n then raise Lsolver.MatrixNotSquare;
          let open Lsolver_impl.Klu in
          sinfo.set_ordering <- klu_set_ordering session;
          sinfo.reinit <- klu_reinit session;
          c_mass_klu session (Matrix.Sparse.sformat smat) m nnz;
          (match sinfo.ordering with None -> ()
                                   | Some o -> c_klu_set_ordering session o)

      | Lsolver_impl.Direct.Superlumt sinfo ->
          if not Sundials_config.superlumt_enabled
            then raise Sundials.NotImplementedBySundialsVersion;
          let smat = Matrix.unwrap mat in
          let m, n = Matrix.Sparse.size smat in
          let nnz, _ = Matrix.Sparse.dims smat in
          if m <> n then raise Lsolver.MatrixNotSquare;
          let open Lsolver_impl.Superlumt in
          sinfo.set_ordering <- superlumt_set_ordering session;
          c_mass_superlumt session m nnz sinfo.num_threads;
          (match sinfo.ordering with None -> ()
                                   | Some o -> c_superlumt_set_ordering session o)

    | Lsolver_impl.Direct.Custom _ ->
        assert false

    let set_mass_callbacks (type m) (type tag)
          (massfn : m mass_fn)
          (solver : (m, 'nd, 'nk, tag) Lsolver_impl.Direct.solver)
          (mat : ('mk, m, 'nd, 'nk) Matrix.t) session =
      let cb = { massfn = massfn; mmat  = (None : m option) } in
      let open Lsolver_impl.Direct in
      begin match solver with
      | Dense ->
          session.mass_callbacks <- DlsDenseMassCallback (cb, Matrix.unwrap mat)
      | LapackDense ->
          session.mass_callbacks <- DlsDenseMassCallback (cb, Matrix.unwrap mat)
      | Band ->
          session.mass_callbacks <- DlsBandMassCallback (cb, Matrix.unwrap mat)
      | LapackBand ->
          session.mass_callbacks <- DlsBandMassCallback (cb, Matrix.unwrap mat)
      | Klu _ ->
          session.mass_callbacks <- SlsKluMassCallback (cb, Matrix.unwrap mat)
      | Superlumt _ ->
          session.mass_callbacks
            <- SlsSuperlumtMassCallback (cb, Matrix.unwrap mat)
      | Custom _ ->
          session.mass_callbacks
            <- DirectCustomMassCallback (cb, Matrix.unwrap mat)
      end;
      session.mass_precfns <- NoMassPrecFns

    (* Sundials >= 3.0.0 *)
    external c_dls_set_mass_linear_solver
      : 'k serial_session
      -> ('m, Nvector_serial.data, 'k) Lsolver_impl.Direct.cptr
        -> ('mk, 'm, Nvector_serial.data, 'k) Matrix.t
        -> bool
        -> unit
      = "c_arkode_dls_set_mass_linear_solver"

    let make { Lsolver_impl.Direct.rawptr; Lsolver_impl.Direct.solver }
             massfn time_dep mat session nv =
      set_mass_callbacks massfn solver mat session;
      if in_compat_mode then make_compat solver mat session
      else c_dls_set_mass_linear_solver session rawptr mat time_dep

    (* Sundials < 3.0.0 *)
    let invalidate_callback session =
      if in_compat_mode then
        match session.mass_callbacks with
        | DlsDenseMassCallback ({ mmat = Some d } as cb, _) ->
            Matrix.Dense.invalidate d;
            cb.mmat <- None
        | DlsBandMassCallback  ({ mmat = Some d } as cb, _) ->
            Matrix.Band.invalidate d;
            cb.mmat <- None
        | SlsKluMassCallback ({ mmat = Some d } as cb, _) ->
            Matrix.Sparse.invalidate d;
            cb.mmat <- None
        | SlsSuperlumtMassCallback ({ mmat = Some d } as cb, _) ->
            Matrix.Sparse.invalidate d;
            cb.mmat <- None
        | _ -> ()

    external get_work_space : 'k serial_session -> int * int
        = "c_arkode_dls_get_mass_work_space"

    let get_work_space s =
      mass_check_direct s;
      get_work_space s

    external c_get_num_mass_setups : 'k serial_session -> int
      = "c_arkode_dls_get_num_mass_setups"

    let get_num_setups s =
      mass_check_direct s;
      if in_compat_mode then raise Sundials.NotImplementedBySundialsVersion else
      c_get_num_mass_setups s

    (* Sundials < 3.0.0 *)
    external c_klu_get_num_mass_evals : 'k serial_session -> int
      = "c_arkode_klu_get_num_mass_evals"

    (* Sundials < 3.0.0 *)
    external c_superlumt_get_num_mass_evals : 'k serial_session -> int
      = "c_arkode_superlumt_get_num_mass_evals"

    external c_get_num_mass_solves : 'k serial_session -> int
      = "c_arkode_dls_get_num_mass_solves"

    let compat_get_num_mass_evals s =
      match s.mass_callbacks with
      | SlsKluMassCallback _ ->
          c_klu_get_num_mass_evals s
      | SlsSuperlumtMassCallback _ -> c_superlumt_get_num_mass_evals s
      | _ -> c_get_num_mass_solves s

    let get_num_solves s =
      mass_check_direct s;
      if in_compat_mode then compat_get_num_mass_evals s else
      c_get_num_mass_solves s

    external c_get_num_mass_mult : 'k serial_session -> int
      = "c_arkode_dls_get_num_mass_mult"

    let get_num_mult s =
      mass_check_direct s;
      if in_compat_mode then raise Sundials.NotImplementedBySundialsVersion else
      c_get_num_mass_mult s

  end (* }}} *)

  module Spils = struct (* {{{ *)
    include MassTypes.Iterative'

    (* Sundials < 3.0.0 *)
    external c_spgmr
      : ('a, 'k) session
        -> int -> Lsolver_impl.Iterative.preconditioning_type -> unit
      = "c_arkode_spils_mass_spgmr"

    (* Sundials < 3.0.0 *)
    external c_spbcgs
      : ('a, 'k) session
        -> int -> Lsolver_impl.Iterative.preconditioning_type -> unit
      = "c_arkode_spils_mass_spbcg"

    (* Sundials < 3.0.0 *)
    external c_sptfqmr
      : ('a, 'k) session
        -> int -> Lsolver_impl.Iterative.preconditioning_type -> unit
      = "c_arkode_spils_mass_sptfqmr"

    (* Sundials < 3.0.0 *)
    external c_spfgmr
      : ('a, 'k) session
        -> int -> Lsolver_impl.Iterative.preconditioning_type -> unit
      = "c_arkode_spils_mass_spfgmr"

    (* Sundials < 3.0.0 *)
    external c_pcg
      : ('a, 'k) session
        -> int -> Lsolver_impl.Iterative.preconditioning_type -> unit
      = "c_arkode_spils_mass_pcg"

    (* Sundials < 3.0.0 *)
    external c_set_gs_type
      : ('a, 'k) session
        -> Lsolver_impl.Iterative.gramschmidt_type -> unit
      = "c_arkode_spils_set_mass_gs_type"

    (* Sundials < 3.0.0 *)
    external c_set_maxl : ('a, 'k) session -> int -> unit
      = "c_arkode_spils_set_mass_maxl"

    (* Sundials < 3.0.0 *)
    external c_set_prec_type
      : ('a, 'k) session -> Lsolver_impl.Iterative.preconditioning_type -> unit
      = "c_arkode_spils_set_mass_prec_type"

    let old_set_gs_type s t =
      mass_check_spils s;
      c_set_gs_type s t

    let old_set_maxl s maxl =
      mass_check_spils s;
      c_set_maxl s maxl

    let old_set_prec_type s t =
      mass_check_spils s;
      c_set_prec_type s t

  external c_set_mass_times : ('a, 'k) session -> bool -> unit
    = "c_arkode_spils_set_mass_times"

    external c_set_preconditioner
      : ('a, 'k) session -> bool -> unit
      = "c_arkode_spils_set_mass_preconditioner"

    external c_spils_set_mass_linear_solver
      : ('a, 'k) session
        -> ('a, 'k) Lsolver_impl.Iterative.cptr
        -> bool
        -> bool
        -> unit
      = "c_arkode_spils_set_mass_linear_solver"

    let init_preconditioner solve setup session nv =
      c_set_preconditioner session (setup <> None);
      session.mass_precfns <- MassPrecFns { prec_solve_fn = solve;
                                            prec_setup_fn = setup }

    let prec_none = Lsolver_impl.Iterative.(PrecNone,
                      fun session nv -> session.mass_precfns <- NoMassPrecFns)

    let prec_left ?setup solve  = Lsolver_impl.Iterative.(PrecLeft,
                                              init_preconditioner solve setup)

    let prec_right ?setup solve = Lsolver_impl.Iterative.(PrecRight,
                                              init_preconditioner solve setup)

    let prec_both ?setup solve  = Lsolver_impl.Iterative.(PrecBoth,
                                              init_preconditioner solve setup)

    let make (type s)
          ({ Lsolver_impl.Iterative.rawptr;
             Lsolver_impl.Iterative.solver;
             Lsolver_impl.Iterative.compat =
               ({ Lsolver_impl.Iterative.maxl;
                  Lsolver_impl.Iterative.gs_type } as compat) })
          ?mass_times_setup mass_times_vec time_dep (prec_type, set_prec)
          session nv =
      if in_compat_mode then begin
        if mass_times_setup <> None then
          raise Sundials.NotImplementedBySundialsVersion;
        let open Lsolver_impl.Iterative in
        (match (solver : ('nd, 'nk, s) solver) with
         | Spgmr ->
             c_spgmr session maxl prec_type;
             (match gs_type with None -> () | Some t -> c_set_gs_type session t);
             compat.set_gs_type <- old_set_gs_type session;
             compat.set_prec_type <- old_set_prec_type session
         | Spfgmr ->
             c_spfgmr session maxl prec_type;
             (match gs_type with None -> () | Some t -> c_set_gs_type session t);
             compat.set_gs_type <- old_set_gs_type session;
             compat.set_prec_type <- old_set_prec_type session
         | Spbcgs ->
             c_spbcgs session maxl prec_type;
             compat.set_maxl <- old_set_maxl session;
             compat.set_prec_type <- old_set_prec_type session
         | Sptfqmr ->
             c_sptfqmr session maxl prec_type;
             compat.set_maxl <- old_set_maxl session;
             compat.set_prec_type <- old_set_prec_type session
         | Pcg ->
             c_pcg session maxl prec_type;
             compat.set_maxl <- old_set_maxl session;
             compat.set_prec_type <- old_set_prec_type session
         | Custom _ -> assert false);
        set_prec session nv;
        session.mass_callbacks <- SpilsMassCallback (mass_times_vec, None)
      end else
        c_spils_set_mass_linear_solver session rawptr
          (mass_times_setup <> None) time_dep;
        Lsolver_impl.Iterative.(c_set_prec_type rawptr solver prec_type);
        set_prec session nv;
        session.mass_callbacks <- SpilsMassCallback
                                    (mass_times_vec, mass_times_setup)

    let set_times s ?mass_times_setup mass_times_vec =
      if in_compat_mode && mass_times_setup <> None then
          raise Sundials.NotImplementedBySundialsVersion;
      match s.mass_callbacks with
      | SpilsMassCallback _ ->
          c_set_mass_times s (mass_times_setup <> None);
          s.mass_callbacks <- SpilsMassCallback
                                (mass_times_vec, mass_times_setup)
      | _ -> raise Sundials.InvalidLinearSolver

    let set_preconditioner s ?setup solve =
      match s.mass_callbacks with
      | SpilsMassCallback _ ->
          c_set_preconditioner s (setup <> None);
          s.mass_precfns <- MassPrecFns { prec_setup_fn = setup;
                                          prec_solve_fn = solve }
      | _ -> raise Sundials.InvalidLinearSolver

    external set_eps_lin            : ('a, 'k) session -> float -> unit
        = "c_arkode_spils_set_mass_eps_lin"

    let set_eps_lin s epsl =
      mass_check_spils s;
      set_eps_lin s epsl

    external get_num_lin_iters      : ('a, 'k) session -> int
        = "c_arkode_spils_get_num_mass_iters"

    let get_num_lin_iters s =
      mass_check_spils s;
      get_num_lin_iters s

    external get_num_conv_fails     : ('a, 'k) session -> int
        = "c_arkode_spils_get_num_mass_conv_fails"

    let get_num_conv_fails s =
      mass_check_spils s;
      get_num_conv_fails s

    external get_num_mtsetup_evals   : ('a, 'k) session -> int
      = "c_arkode_spils_get_num_mtsetup_evals"

    let get_num_mtsetup_evals s =
      mass_check_spils s;
      get_num_mtsetup_evals s

    external get_num_mtimes_evals   : ('a, 'k) session -> int
        = "c_arkode_spils_get_num_mtimes_evals"

    let get_num_mtimes_evals s =
      mass_check_spils s;
      get_num_mtimes_evals s

    external get_work_space         : ('a, 'k) session -> int * int
        = "c_arkode_spils_get_mass_work_space"

    let get_work_space s =
      mass_check_spils s;
      get_work_space s

    external get_num_prec_evals     : ('a, 'k) session -> int
        = "c_arkode_spils_get_num_mass_prec_evals"

    let get_num_prec_evals s =
      mass_check_spils s;
      get_num_prec_evals s

    external get_num_prec_solves    : ('a, 'k) session -> int
        = "c_arkode_spils_get_num_mass_prec_solves"

    let get_num_prec_solves s =
      mass_check_spils s;
      get_num_prec_solves s
  end (* }}} *)

  module Alternate = struct (* {{{ *)
    include MassTypes.Alternate'

    external c_set_alternate
      : ('data, 'kind) session -> bool -> bool -> unit
      = "c_arkode_set_mass_alternate"

    let make f s nv =
      let { minit; msetup; msolve } as cb = f s nv in
      c_set_alternate s (minit <> None) (msetup <> None);
      s.mass_precfns <- NoMassPrecFns;
      s.mass_callbacks <- AlternateMassCallback cb
  end (* }}} *)
end (* }}} *)


external sv_tolerances  : ('a, 'k) session -> float -> ('a, 'k) nvector -> unit
    = "c_arkode_sv_tolerances"
external ss_tolerances  : ('a, 'k) session -> float -> float -> unit
    = "c_arkode_ss_tolerances"
external wf_tolerances  : ('a, 'k) session -> unit
    = "c_arkode_wf_tolerances"

type ('a, 'k) tolerance =
  | SStolerances of float * float
  | SVtolerances of float * ('a, 'k) nvector
  | WFtolerances of 'a error_weight_fun

let default_tolerances = SStolerances (1.0e-4, 1.0e-9)

let set_tolerances s tol =
  match tol with
  | SStolerances (rel, abs) -> (s.errw <- dummy_errw; ss_tolerances s rel abs)
  | SVtolerances (rel, abs) -> (if Sundials_config.safe then s.checkvec abs;
                                s.errw <- dummy_errw; sv_tolerances s rel abs)
  | WFtolerances ferrw -> (s.errw <- ferrw; wf_tolerances s)

external ress_tolerance  : ('a, 'k) session -> float -> unit
    = "c_arkode_ress_tolerance"
external resv_tolerance  : ('a, 'k) session -> ('a, 'k) nvector -> unit
    = "c_arkode_resv_tolerance"
external resf_tolerance  : ('a, 'k) session -> unit
    = "c_arkode_resf_tolerance"

type ('data, 'kind) res_tolerance =
  | ResStolerance of float
  | ResVtolerance of ('data, 'kind) Nvector.t
  | ResFtolerance of 'data res_weight_fun

let set_res_tolerance s tol =
  match tol with
  | ResStolerance abs -> (s.resw <- dummy_resw;
                          s.uses_resv <- false;
                          ress_tolerance s abs)
  | ResVtolerance abs -> (if Sundials_config.safe then s.checkvec abs;
                          s.uses_resv <- true;
                          s.resw <- dummy_resw;
                          resv_tolerance s abs)
  | ResFtolerance fresw -> (s.resw <- fresw;
                            s.uses_resv <- false;
                            resf_tolerance s)

external set_fixed_point : ('a, 'k) session -> int -> unit
  = "c_arkode_set_fixed_point"

external set_newton : ('a, 'k) session -> unit
  = "c_arkode_set_newton"

external set_linear : ('a, 'k) session -> bool -> unit
  = "c_arkode_set_linear"

external set_nonlinear : ('a, 'k) session -> unit
  = "c_arkode_set_nonlinear"

external c_set_order : ('a, 'k) session -> int -> unit
  = "c_arkode_set_order"

external c_session_finalize : ('a, 'kind) session -> unit
    = "c_arkode_session_finalize"

let session_finalize s =
  Dls.invalidate_callback s;
  c_session_finalize s

external c_init :
  ('a, 'k) session Weak.t
  -> bool             (* f_i given *)
  -> bool             (* f_e given *)
  -> ('a, 'k) nvector (* y_0 *)
  -> float            (* t_0 *)
  -> (arkode_mem * c_weak_ref)
  = "c_arkode_init"

let init prob tol ?restol ?order ?mass ?(roots=no_roots) t0 y0 =
  let (nroots, roots) = roots in
  let checkvec = Nvector.check y0 in
  if Sundials_config.safe && nroots < 0 then
    raise (Invalid_argument "number of root functions is negative");
  let weakref = Weak.create 1 in
  let problem, fi, fe, iter, lin =
    match prob with
    | Implicit (fi,i,l) -> ImplicitOnly, Some fi, None,    Some i, Some l
    | Explicit fe       -> ExplicitOnly, None,    Some fe, None,   None
    | ImEx { implicit=(fi, i, l); explicit=fe }
                        -> ImplicitAndExplicit, Some fi, Some fe, Some i, Some l
  in
  let arkode_mem, backref = c_init weakref (fi <> None) (fe <> None) y0 t0 in
  (* arkode_mem and backref have to be immediately captured in a session and
     associated with the finalizer before we do anything else.  *)
  let session = {
          arkode       = arkode_mem;
          backref      = backref;
          nroots       = nroots;
          checkvec     = checkvec;
          uses_resv    = false;

          exn_temp     = None;

          problem      = problem;
          irhsfn       = (match fi with Some f -> f | None -> dummy_irhsfn);
          erhsfn       = (match fe with Some f -> f | None -> dummy_erhsfn);

          rootsfn      = roots;
          errh         = dummy_errh;
          errw         = dummy_errw;
          resw         = dummy_resw;

          adaptfn      = dummy_adaptfn;
          stabfn       = dummy_stabfn;
          resizefn     = dummy_resizefn;
          poststepfn   = dummy_poststepfn;

          linsolver      = None;
          ls_solver      = Lsolver_impl.NoSolver;
          ls_callbacks   = NoCallbacks;
          ls_precfns     = NoPrecFns;
          mass_solver    = Lsolver_impl.NoSolver;
          mass_callbacks = NoMassCallbacks;
          mass_precfns   = NoMassPrecFns;
        } in
  Gc.finalise session_finalize session;
  Weak.set weakref 0 (Some session);
  (* Now the session is safe to use.  If any of the following fails and raises
     an exception, the GC will take care of freeing arkode_mem and backref.  *)
  if nroots > 0 then
    c_root_init session nroots;
  (match iter with
   | None -> ()
   | Some (Newton ls) -> (session.linsolver <- Some ls;
                          ls session y0)
   | Some (FixedPoint fpm) -> set_fixed_point session fpm);
  (match lin with
   | Some (Linear timedepend) -> set_linear session timedepend
   | _ -> ());
  set_tolerances session tol;
  (match restol with Some rtol -> set_res_tolerance session rtol | None -> ());
  (match order with Some o -> c_set_order session o | None -> ());
  (match mass with Some msolver -> msolver session y0 | None -> ());
  session

let get_num_roots { nroots } = nroots

external c_reinit
    : ('a, 'k) session -> float -> ('a, 'k) nvector -> unit
    = "c_arkode_reinit"

let reinit session ?problem ?order ?roots t0 y0 =
  if Sundials_config.safe then session.checkvec y0;
  Dls.invalidate_callback session;
  (match problem with
   | None -> ()
   | Some (Implicit (fi, _, _)) ->
       (session.problem <- ImplicitOnly;
        session.irhsfn <- fi;
        session.erhsfn <- dummy_erhsfn)
   | Some (Explicit fe) ->
       (session.problem <- ExplicitOnly;
        session.irhsfn <- dummy_irhsfn;
        session.erhsfn <- fe)
   | Some (ImEx { implicit=(fi, _, _); explicit=fe }) ->
       (session.problem <- ImplicitAndExplicit;
        session.irhsfn <- fi;
        session.erhsfn <- fe));
  c_reinit session t0 y0;
  (match order with Some o -> c_set_order session o | None -> ());
  (match roots with
   | None -> ()
   | Some roots -> root_init session roots)

external c_resize
    : ('a, 'k) session -> bool -> float -> float -> ('a, 'k) nvector -> unit
    = "c_arkode_resize"

let resize session ?resize_nvec ?linsolv tol ?restol hscale ynew t0 =
  session.checkvec <- Nvector.check ynew;
  (match linsolv with None -> () | ls -> session.linsolver <- ls);
  (match resize_nvec with None -> () | Some f -> session.resizefn <- f);
  c_resize session (resize_nvec <> None) hscale t0 ynew;
  session.resizefn <- dummy_resizefn;
  (match Sundials.sundials_version with
   | 2,6,1 | 2,6,2 -> () (* avoid a segmentation fault in earlier versions *)
   | _ -> set_tolerances session tol);
  (match restol with
   | Some rt -> set_res_tolerance session rt
   | _ when session.uses_resv -> set_res_tolerance session (ResStolerance 1.e-9)
   | _ -> ());
  (match session.linsolver with Some ls -> ls session ynew | None -> ())

external get_root_info  : ('a, 'k) session -> Sundials.Roots.t -> unit
    = "c_arkode_get_root_info"

type solver_result =
  | Success             (** ARK_SUCCESS *)
  | RootsFound          (** ARK_ROOT_RETURN *)
  | StopTimeReached     (** ARK_TSTOP_RETURN *)

external c_solve_normal : ('a, 'k) session -> float -> ('a, 'k) nvector
                              -> float * solver_result
    = "c_arkode_solve_normal"

let solve_normal s t y =
  if Sundials_config.safe then s.checkvec y;
  c_solve_normal s t y

external c_solve_one_step : ('a, 'k) session -> float -> ('a, 'k) nvector
                              -> float * solver_result
    = "c_arkode_solve_one_step"

let solve_one_step s t y =
  if Sundials_config.safe then s.checkvec y;
  c_solve_one_step s t y

external c_get_dky
    : ('a, 'k) session -> float -> int -> ('a, 'k) nvector -> unit
    = "c_arkode_get_dky"

let get_dky s y =
  if Sundials_config.safe then s.checkvec y;
  fun t k -> c_get_dky s t k y

external get_integrator_stats : ('a, 'k) session -> integrator_stats
    = "c_arkode_get_integrator_stats"

external get_work_space         : ('a, 'k) session -> int * int
    = "c_arkode_get_work_space"

external get_num_steps          : ('a, 'k) session -> int
    = "c_arkode_get_num_steps"

external get_num_exp_steps      : ('d, 'k) session -> int
    = "c_arkode_get_num_exp_steps"

external get_num_acc_steps      : ('d, 'k) session -> int
    = "c_arkode_get_num_acc_steps"

external get_num_step_attempts  : ('d, 'k) session -> int
    = "c_arkode_get_num_step_attempts"

external get_num_rhs_evals      : ('a, 'k) session -> int * int
    = "c_arkode_get_num_rhs_evals"

external get_num_lin_solv_setups : ('a, 'k) session -> int
    = "c_arkode_get_num_lin_solv_setups"

external get_num_mass_solves     : ('d, 'k) session -> int
    = "c_arkode_get_num_mass_solves"

external get_num_err_test_fails : ('a, 'k) session -> int
    = "c_arkode_get_num_err_test_fails"

external get_actual_init_step   : ('a, 'k) session -> float
    = "c_arkode_get_actual_init_step"

external get_last_step          : ('a, 'k) session -> float
    = "c_arkode_get_last_step"

external get_current_step       : ('a, 'k) session -> float
    = "c_arkode_get_current_step"

external get_current_time       : ('a, 'k) session -> float
    = "c_arkode_get_current_time"

let print_integrator_stats s oc =
  let stats = get_integrator_stats s
  in
    Printf.fprintf oc "num_steps = %d\n"           stats.num_steps;
    Printf.fprintf oc "exp_steps = %d\n"           stats.exp_steps;
    Printf.fprintf oc "acc_steps = %d\n"           stats.acc_steps;
    Printf.fprintf oc "step_attempts = %d\n"       stats.step_attempts;
    Printf.fprintf oc "num_nfe_evals = %d\n"       stats.num_nfe_evals;
    Printf.fprintf oc "num_nfi_evals = %d\n"       stats.num_nfi_evals;
    Printf.fprintf oc "num_lin_solv_setups = %d\n" stats.num_lin_solv_setups;
    Printf.fprintf oc "num_err_test_fails = %d\n"  stats.num_err_test_fails;
    Printf.fprintf oc "actual_init_step = %e\n"    stats.actual_init_step;
    Printf.fprintf oc "last_step = %e\n"           stats.last_step;
    Printf.fprintf oc "current_step = %e\n"        stats.current_step;
    Printf.fprintf oc "current_time = %e\n"        stats.current_time

external set_diagnostics : ('a, 'k) session -> Sundials.Logfile.t -> unit
    = "c_arkode_set_diagnostics"

external clear_diagnostics : ('a, 'k) session -> unit
    = "c_arkode_clear_diagnostics"

external set_error_file : ('a, 'k) session -> Sundials.Logfile.t -> unit
    = "c_arkode_set_error_file"

external c_set_err_handler_fn  : ('a, 'k) session -> unit
    = "c_arkode_set_err_handler_fn"

let set_err_handler_fn s ferrh =
  s.errh <- ferrh;
  c_set_err_handler_fn s

external clear_err_handler_fn  : ('a, 'k) session -> unit
    = "c_arkode_clear_err_handler_fn"

let clear_err_handler_fn s =
  s.errh <- dummy_errh;
  clear_err_handler_fn s

external c_set_imex             : ('a, 'k) session -> unit
    = "c_arkode_set_imex"

external c_set_explicit         : ('a, 'k) session -> unit
    = "c_arkode_set_explicit"

external c_set_implicit         : ('a, 'k) session -> unit
    = "c_arkode_set_implicit"

let set_imex s =
  (if s.irhsfn == dummy_irhsfn || s.erhsfn == dummy_erhsfn then raise IllInput);
  s.problem <- ImplicitAndExplicit;
  c_set_imex s

let set_explicit s =
  (if s.erhsfn == dummy_erhsfn then raise IllInput);
  s.problem <- ExplicitOnly;
  c_set_explicit s

let set_implicit s =
  (if s.irhsfn == dummy_irhsfn then raise IllInput);
  s.problem <- ImplicitOnly;
  c_set_implicit s

type rk_method = {
    stages : int;
    global_order : int;
    global_embedded_order : int;
  }

type rk_timescoefs = {
    stage_times  : RealArray.t;
    coefficients : RealArray.t;
    bembed       : RealArray.t option;
  }

external c_set_ark_tables
  : ('d, 'k) session -> rk_method
    -> RealArray.t -> RealArray.t
    -> rk_timescoefs * rk_timescoefs
    -> unit
    = "c_arkode_set_ark_tables"

external c_set_erk_table
  : ('d, 'k) session -> rk_method -> RealArray.t -> rk_timescoefs -> unit
    = "c_arkode_set_erk_table"

external c_set_irk_table
  : ('d, 'k) session -> rk_method -> RealArray.t -> rk_timescoefs -> unit
    = "c_arkode_set_irk_table"

let set_ark_tables s rkm ai ae tci tce =
  (if s.irhsfn == dummy_irhsfn || s.erhsfn == dummy_erhsfn then raise IllInput);
  (match Sundials.sundials_version with
   | 2,5,_ | 2,6,_ -> if tci.bembed = None || tce.bembed = None
                      then raise Sundials.NotImplementedBySundialsVersion
   | _ -> ());
  c_set_ark_tables s rkm ai ae (tci, tce)

let set_erk_table s rkm ae tc =
  (if s.erhsfn == dummy_erhsfn then raise IllInput);
  (match Sundials.sundials_version with
   | 2,5,_ | 2,6,_ -> if tc.bembed = None
                      then raise Sundials.NotImplementedBySundialsVersion
   | _ -> ());
  c_set_erk_table s rkm ae tc

let set_irk_table s rkm ai tc =
  (if s.irhsfn == dummy_irhsfn then raise IllInput);
  (match Sundials.sundials_version with
   | 2,5,_ | 2,6,_ -> if tc.bembed = None
                      then raise Sundials.NotImplementedBySundialsVersion
   | _ -> ());
  c_set_irk_table s rkm ai tc

type erk_table =
  | HeunEuler_2_1_2
  | BogackiShampine_4_2_3
  | ARK_4_2_3_Explicit
  | Zonneveld_5_3_4
  | ARK_6_3_4_Explicit
  | SayfyAburub_6_3_4
  | CashKarp_6_4_5
  | Fehlberg_6_4_5
  | DormandPrince_7_4_5
  | ARK_8_4_5_Explicit
  | Verner_8_5_6
  | Fehlberg_13_7_8

type irk_table =
  | SDIRK_2_1_2
  | Billington_3_2_3
  | TRBDF2_3_2_3
  | Kvaerno_4_2_3
  | ARK_4_2_3_Implicit
  | Cash_5_2_4
  | Cash_5_3_4
  | SDIRK_5_3_4
  | Kvaerno_5_3_4
  | ARK_6_3_4_Implicit
  | Kvaerno_7_4_5
  | ARK_8_4_5_Implicit

type ark_table =
  | ARK_4_2_3
  | ARK_6_3_4
  | ARK_8_4_5

let int_of_erk_table v =
  match Sundials.sundials_version with
  | 2,5,_ | 2,6,_ ->
    (match v with
     | HeunEuler_2_1_2       -> 0
     | BogackiShampine_4_2_3 -> 1
     | ARK_4_2_3_Explicit    -> 2
     | Zonneveld_5_3_4       -> 3
     | ARK_6_3_4_Explicit    -> 4
     | SayfyAburub_6_3_4     -> 5
     | CashKarp_6_4_5        -> 6
     | Fehlberg_6_4_5        -> 7
     | DormandPrince_7_4_5   -> 8
     | ARK_8_4_5_Explicit    -> 9
     | Verner_8_5_6          -> 10
     | Fehlberg_13_7_8       -> raise Sundials.NotImplementedBySundialsVersion)
  | _ ->
    (match v with
     | HeunEuler_2_1_2       -> 0
     | BogackiShampine_4_2_3 -> 1
     | ARK_4_2_3_Explicit    -> 2
     | Zonneveld_5_3_4       -> 3
     | ARK_6_3_4_Explicit    -> 4
     | SayfyAburub_6_3_4     -> 5
     | CashKarp_6_4_5        -> 6
     | Fehlberg_6_4_5        -> 7
     | DormandPrince_7_4_5   -> 8
     | ARK_8_4_5_Explicit    -> 9
     | Verner_8_5_6          -> 10
     | Fehlberg_13_7_8       -> 11)

let int_of_irk_table v =
  match Sundials.sundials_version with
  | 2,5,_ | 2,6,_ ->
    (match v with
     | SDIRK_2_1_2        -> 11
     | Billington_3_2_3   -> 12
     | TRBDF2_3_2_3       -> 13
     | Kvaerno_4_2_3      -> 14
     | ARK_4_2_3_Implicit -> 15
     | Cash_5_2_4         -> 16
     | Cash_5_3_4         -> 17
     | SDIRK_5_3_4        -> 18
     | Kvaerno_5_3_4      -> 19
     | ARK_6_3_4_Implicit -> 20
     | Kvaerno_7_4_5      -> 21
     | ARK_8_4_5_Implicit -> 22)
  | _ ->
    (match v with
     | SDIRK_2_1_2        -> 12
     | Billington_3_2_3   -> 13
     | TRBDF2_3_2_3       -> 14
     | Kvaerno_4_2_3      -> 15
     | ARK_4_2_3_Implicit -> 16
     | Cash_5_2_4         -> 17
     | Cash_5_3_4         -> 18
     | SDIRK_5_3_4        -> 19
     | Kvaerno_5_3_4      -> 20
     | ARK_6_3_4_Implicit -> 21
     | Kvaerno_7_4_5      -> 22
     | ARK_8_4_5_Implicit -> 23)

let ints_of_ark_table v =
  match v with
  | ARK_4_2_3 -> (int_of_irk_table ARK_4_2_3_Implicit,
                  int_of_erk_table ARK_4_2_3_Explicit)
  | ARK_6_3_4 -> (int_of_irk_table ARK_6_3_4_Implicit,
                  int_of_erk_table ARK_6_3_4_Explicit)
  | ARK_8_4_5 -> (int_of_irk_table ARK_8_4_5_Implicit,
                  int_of_erk_table ARK_8_4_5_Explicit)

external c_set_erk_table_num      : ('d, 'k) session -> int -> unit
    = "c_arkode_set_erk_table_num"
external c_set_irk_table_num      : ('d, 'k) session -> int -> unit
    = "c_arkode_set_irk_table_num"
external c_set_ark_table_num      : ('d, 'k) session -> int * int -> unit
    = "c_arkode_set_ark_table_num"

let set_erk_table_num s v = c_set_erk_table_num s (int_of_erk_table v)
let set_irk_table_num s v = c_set_irk_table_num s (int_of_irk_table v)
let set_ark_table_num s v = c_set_ark_table_num s (ints_of_ark_table v)

type adaptivity_params = {
    ks : (float * float * float) option;
    method_order : bool;
  }

type 'd adaptivity_method =
  | PIDcontroller of adaptivity_params
  | PIcontroller of adaptivity_params
  | Icontroller of adaptivity_params
  | ExplicitGustafsson of adaptivity_params
  | ImplicitGustafsson of adaptivity_params
  | ImExGustafsson of adaptivity_params
  | AdaptivityFn of 'd adaptivity_fn

external c_set_adaptivity_method
    : ('d, 'k) session -> 'd adaptivity_method -> unit
    = "c_arkode_set_adaptivity_method"

let set_adaptivity_method s am =
  (match am with
   | AdaptivityFn fn -> s.adaptfn <- fn
   | _ -> s.adaptfn <- dummy_adaptfn);
  c_set_adaptivity_method s am

external c_set_stability_fn : ('d, 'k) session -> bool -> unit
    = "c_arkode_set_stability_fn"

let set_stability_fn s f =
  s.stabfn <- f;
  c_set_stability_fn s true

let clear_stability_fn s =
  s.stabfn <- dummy_stabfn;
  c_set_stability_fn s false

type predictor_method =
  | TrivialPredictor
  | MaximumOrderPredictor
  | VariableOrderPredictor
  | CutoffOrderPredictor
  | BootstrapPredictor

external set_predictor_method : ('d, 'k) session -> predictor_method -> unit
    = "c_arkode_set_predictor_method"

external set_defaults           : ('a, 'k) session -> unit
    = "c_arkode_set_defaults"
external set_dense_order        : ('a, 'k) session -> int -> unit
    = "c_arkode_set_dense_order"
external set_max_num_steps      : ('a, 'k) session -> int -> unit
    = "c_arkode_set_max_num_steps"
external set_max_hnil_warns     : ('a, 'k) session -> int -> unit
    = "c_arkode_set_max_hnil_warns"
external set_init_step          : ('a, 'k) session -> float -> unit
    = "c_arkode_set_init_step"
external set_fixed_step         : ('a, 'k) session -> float option -> unit
    = "c_arkode_set_fixed_step"
external set_min_step           : ('a, 'k) session -> float -> unit
    = "c_arkode_set_min_step"
external set_max_step           : ('a, 'k) session -> float -> unit
    = "c_arkode_set_max_step"
external set_stop_time          : ('a, 'k) session -> float -> unit
    = "c_arkode_set_stop_time"
external set_optimal_params     : ('a, 'k) session -> unit
    = "c_arkode_set_optimal_params"
external set_max_err_test_fails : ('a, 'k) session -> int -> unit
    = "c_arkode_set_max_err_test_fails"
external set_max_nonlin_iters   : ('a, 'k) session -> int -> unit
    = "c_arkode_set_max_nonlin_iters"
external set_max_conv_fails     : ('a, 'k) session -> int -> unit
    = "c_arkode_set_max_conv_fails"
external set_nonlin_conv_coef   : ('a, 'k) session -> float -> unit
    = "c_arkode_set_nonlin_conv_coef"
external set_nonlin_crdown      : ('a, 'k) session -> float -> unit
    = "c_arkode_set_nonlin_crdown"
external set_nonlin_rdiv        : ('a, 'k) session -> float -> unit
    = "c_arkode_set_nonlin_rdiv"
external set_delta_gamma_max    : ('a, 'k) session -> float -> unit
    = "c_arkode_set_delta_gamma_max"
external set_max_steps_between_lset : ('a, 'k) session -> int -> unit
    = "c_arkode_set_max_steps_between_lset"
external set_cfl_fraction       : ('a, 'k) session -> float -> unit
    = "c_arkode_set_cfl_fraction"
external set_error_bias         : ('a, 'k) session -> float -> unit
    = "c_arkode_set_error_bias"
external set_fixed_step_bounds  : ('a, 'k) session -> float -> float -> unit
    = "c_arkode_set_fixed_step_bounds"
external set_max_cfail_growth   : ('a, 'k) session -> float -> unit
    = "c_arkode_set_max_cfail_growth"
external set_max_efail_growth   : ('a, 'k) session -> float -> unit
    = "c_arkode_set_max_efail_growth"
external set_max_first_growth   : ('a, 'k) session -> float -> unit
    = "c_arkode_set_max_first_growth"
external set_max_growth         : ('a, 'k) session -> float -> unit
    = "c_arkode_set_max_growth"
external set_safety_factor      : ('a, 'k) session -> float -> unit
    = "c_arkode_set_safety_factor"
external set_small_num_efails   : ('a, 'k) session -> float -> unit
    = "c_arkode_set_small_num_efails"

external c_set_postprocess_step_fn : ('a, 'k) session -> bool -> unit
    = "c_arkode_set_postprocess_step_fn"

let set_postprocess_step_fn s fn =
  s.poststepfn <- fn;
  c_set_postprocess_step_fn s true

let clear_postprocess_step_fn s =
  s.poststepfn <- dummy_poststepfn;
  c_set_postprocess_step_fn s false

external c_set_root_direction   : ('a, 'k) session -> Sundials.RootDirs.t -> unit
    = "c_arkode_set_root_direction"

let set_root_direction s rda =
  c_set_root_direction s (Sundials.RootDirs.copy (get_num_roots s) rda)

let set_all_root_directions s rd =
  c_set_root_direction s (Sundials.RootDirs.make (get_num_roots s) rd)

external set_no_inactive_root_warn      : ('a, 'k) session -> unit
    = "c_arkode_set_no_inactive_root_warn"

external get_current_butcher_tables
  : ('d, 'k) session
    -> rk_method * RealArray.t * RealArray.t * rk_timescoefs * rk_timescoefs
    = "c_arkode_get_current_butcher_tables"

external get_tol_scale_factor           : ('a, 'k) session -> float
    = "c_arkode_get_tol_scale_factor"

external c_get_err_weights : ('a, 'k) session -> ('a, 'k) nvector -> unit
    = "c_arkode_get_err_weights"

let get_err_weights s ew =
  if Sundials_config.safe then s.checkvec ew;
  c_get_err_weights s ew

external c_get_est_local_errors : ('a, 'k) session -> ('a, 'k) nvector -> unit
    = "c_arkode_get_est_local_errors"

let get_est_local_errors s ew =
  if Sundials_config.safe then s.checkvec ew;
  c_get_est_local_errors s ew

external get_num_nonlin_solv_iters      : ('a, 'k) session -> int
    = "c_arkode_get_num_nonlin_solv_iters"

external get_num_nonlin_solv_conv_fails : ('a, 'k) session -> int
    = "c_arkode_get_num_nonlin_solv_conv_fails"

external get_nonlin_solv_stats          : ('a, 'k) session -> int * int
    = "c_arkode_get_nonlin_solv_stats"

external get_num_g_evals                : ('a, 'k) session -> int
    = "c_arkode_get_num_g_evals"


(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "c_arkode_init_module"

let _ =
  c_init_module
    (* Exceptions must be listed in the same order as
       arkode_exn_index.  *)
    [|IllInput;
      TooClose;
      TooMuchWork;
      TooMuchAccuracy;
      ErrFailure;
      ConvergenceFailure;
      LinearInitFailure;
      LinearSetupFailure;
      LinearSolveFailure;
      MassInitFailure;
      MassSetupFailure;
      MassSolveFailure;
      MassMultFailure;
      RhsFuncFailure;
      FirstRhsFuncFailure;
      RepeatedRhsFuncFailure;
      UnrecoverableRhsFuncFailure;
      RootFuncFailure;
      PostprocStepFailure;
      BadK;
      BadT;
    |]
