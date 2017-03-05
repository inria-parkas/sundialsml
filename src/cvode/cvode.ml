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

    let solver session nv =
      c_cvode_diag session;
      session.ls_precfns <- NoPrecFns;
      session.ls_callbacks <- DiagNoCallbacks

    external get_work_space       : ('a, 'k) session -> int * int
        = "c_cvode_diag_get_work_space"

    let get_work_space s =
      ls_check_diag s;
      get_work_space s

    external get_num_rhs_evals    : ('a, 'k) session -> int
        = "c_cvode_diag_get_num_rhs_evals"

    let get_num_rhs_evals s =
      ls_check_diag s;
      get_num_rhs_evals s
  end

module Dls =
  struct
    include DlsTypes

    external c_dls_dense : 'k serial_session -> int -> bool -> unit
      = "c_cvode_dls_dense"

    external c_dls_lapack_dense : 'k serial_session -> int -> bool -> unit
      = "c_cvode_dls_lapack_dense"

    external c_dls_band : ('k serial_session * int) -> int -> int -> bool -> unit
      = "c_cvode_dls_band"

    external c_dls_lapack_band
      : ('k serial_session * int) -> int -> int -> bool -> unit
      = "c_cvode_dls_lapack_band"

    external set_dense_jac_fn : 'k serial_session -> unit
        = "c_cvode_dls_set_dense_jac_fn"

    let dense ?jac () session nv =
      let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
      (session.ls_callbacks <-
        match jac with
        | None   -> DlsDenseCallback no_dense_callback
        | Some f -> DlsDenseCallback { jacfn = f; dmat = None });
      session.ls_precfns <- NoPrecFns;
      c_dls_dense session neqs (jac <> None)

    let lapack_dense ?jac () session nv =
      let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
      (session.ls_callbacks <-
        match jac with
        | None   -> DlsDenseCallback no_dense_callback
        | Some f -> DlsDenseCallback { jacfn = f; dmat = None });
      session.ls_precfns <- NoPrecFns;
      c_dls_lapack_dense session neqs (jac <> None)

    let band ?jac p session nv =
      let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
      (session.ls_callbacks <-
        match jac with
        | None   -> DlsBandCallback no_band_callback
        | Some f -> DlsBandCallback { bjacfn = f; bmat = None });
      session.ls_precfns <- NoPrecFns;
      c_dls_band (session, neqs) p.mupper p.mlower (jac <> None)

    let lapack_band ?jac p session nv =
      let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
      (session.ls_callbacks <-
        match jac with
        | None   -> DlsBandCallback no_band_callback
        | Some f -> DlsBandCallback { bjacfn = f; bmat = None });
      session.ls_precfns <- NoPrecFns;
      c_dls_lapack_band (session, neqs) p.mupper p.mlower (jac <> None)

    let invalidate_callback session =
      match session.ls_callbacks with
      | DlsDenseCallback ({ dmat = Some d } as cb) ->
          Dls.DenseMatrix.invalidate d;
          cb.dmat <- None
      | DlsBandCallback  ({ bmat = Some d } as cb) ->
          Dls.BandMatrix.invalidate d;
          cb.bmat <- None
      | SlsKluCallback ({ SlsTypes.smat = Some d } as cb) ->
          Sls_impl.invalidate d;
          cb.SlsTypes.smat <- None
      | SlsSuperlumtCallback ({ SlsTypes.smat = Some d } as cb) ->
          Sls_impl.invalidate d;
          cb.SlsTypes.smat <- None
      | _ -> ()

    let set_dense_jac_fn s fjacfn =
      match s.ls_callbacks with
      | DlsDenseCallback _ ->
          invalidate_callback s;
          s.ls_callbacks <- DlsDenseCallback { jacfn = fjacfn; dmat = None };
          set_dense_jac_fn s
      | _ -> raise Sundials.InvalidLinearSolver

    external clear_dense_jac_fn : 'k serial_session -> unit
        = "c_cvode_dls_clear_dense_jac_fn"

    let clear_dense_jac_fn s =
      match s.ls_callbacks with
      | DlsDenseCallback _ ->
          invalidate_callback s;
          s.ls_callbacks <- DlsDenseCallback no_dense_callback;
          clear_dense_jac_fn s
      | _ -> raise Sundials.InvalidLinearSolver

    external set_band_jac_fn : 'k serial_session -> unit
        = "c_cvode_dls_set_band_jac_fn"

    let set_band_jac_fn s f =
      match s.ls_callbacks with
      | DlsBandCallback _ ->
          invalidate_callback s;
          s.ls_callbacks <- DlsBandCallback { bjacfn = f; bmat = None };
          set_band_jac_fn s
      | _ -> raise Sundials.InvalidLinearSolver

    external clear_band_jac_fn : 'k serial_session -> unit
        = "c_cvode_dls_clear_band_jac_fn"

    let clear_band_jac_fn s =
      match s.ls_callbacks with
      | DlsBandCallback _ ->
          invalidate_callback s;
          s.ls_callbacks <- DlsBandCallback no_band_callback;
          clear_band_jac_fn s
      | _ -> raise Sundials.InvalidLinearSolver

    external get_work_space : 'k serial_session -> int * int
        = "c_cvode_dls_get_work_space"

    let get_work_space s =
      ls_check_dls s;
      get_work_space s

    external get_num_jac_evals : 'k serial_session -> int
        = "c_cvode_dls_get_num_jac_evals"

    let get_num_jac_evals s =
      ls_check_dls s;
      get_num_jac_evals s

    external get_num_rhs_evals : 'k serial_session -> int
        = "c_cvode_dls_get_num_rhs_evals"

    let get_num_rhs_evals s =
      ls_check_dls s;
      get_num_rhs_evals s
  end

module Sls = struct
  include SlsTypes

  module Klu = struct
    (* Must correspond with cvode_klu_ordering_tag *)
    type ordering =
         Amd
       | ColAmd
       | Natural

    external c_klu : 'k serial_session -> Sls_impl.sformat -> int -> int -> unit
      = "c_cvode_klu_init"

    let solver fmt f nnz session nv =
      if not Sundials_config.klu_enabled
        then raise Sundials.NotImplementedBySundialsVersion;
      let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
      session.ls_callbacks <- SlsKluCallback { jacfn = f; smat = None };
      session.ls_precfns <- NoPrecFns;
      c_klu session fmt neqs nnz

    (* We force the type argument here to avoid propagating it to the
       session type; which is unnecessary and needlessy complicated
       for users. *)
    let solver_csc (f : Sls.SparseMatrix.csc sparse_jac_fn)
      = solver Sls_impl.CSC_MAT (Obj.magic f : unit sparse_jac_fn)

    let solver_csr (f : Sls.SparseMatrix.csr sparse_jac_fn)
      = match Sundials.sundials_version with
        | 2,5,_ | 2,6,_ -> raise Sundials.NotImplementedBySundialsVersion
        | _ -> solver Sls_impl.CSR_MAT (Obj.magic f : unit sparse_jac_fn)

    external c_set_ordering : 'k serial_session -> ordering -> unit
      = "c_cvode_klu_set_ordering"

    let set_ordering session ordering =
      ls_check_klu session;
      c_set_ordering session ordering

    external c_reinit : 'k serial_session -> int -> int -> bool -> unit
      = "c_cvode_klu_reinit"

    let reinit session n nnz realloc =
      ls_check_klu session;
      c_reinit session n nnz realloc

    external c_get_num_jac_evals : 'k serial_session -> int
      = "c_cvode_klu_get_num_jac_evals"

    let get_num_jac_evals session =
      ls_check_klu session;
      c_get_num_jac_evals session
  end
    
  module Superlumt = struct
    (* Must correspond with cvode_superlumt_ordering_tag *)
    type ordering =
         Natural
       | MinDegreeProd
       | MinDegreeSum
       | ColAmd

    external c_superlumt : 'k serial_session -> int -> int -> int -> unit
      = "c_cvode_superlumt_init"

    let solver sformat f ~nnz ~nthreads session nv =
      if not Sundials_config.superlumt_enabled
        then raise Sundials.NotImplementedBySundialsVersion;
      let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
      session.ls_callbacks <- SlsSuperlumtCallback { jacfn = f; smat = None };
      session.ls_precfns <- NoPrecFns;
      c_superlumt session neqs nnz nthreads

    (* We force the type argument here to avoid propagating it to the
       session type; which is unnecessary and needlessy complicated
       for users. *)
    let solver_csc (f : Sls.SparseMatrix.csc sparse_jac_fn)
      = solver Sls_impl.CSC_MAT (Obj.magic f : unit sparse_jac_fn)

    external c_set_ordering : 'k serial_session -> ordering -> unit
      = "c_cvode_superlumt_set_ordering"

    let set_ordering session ordering =
      ls_check_superlumt session;
      c_set_ordering session ordering

    external c_get_num_jac_evals : 'k serial_session -> int
      = "c_cvode_superlumt_get_num_jac_evals"

    let get_num_jac_evals session =
      ls_check_superlumt session;
      c_get_num_jac_evals session
  end
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

    let init_preconditioner solve setup session nv =
      c_set_preconditioner session (setup <> None);
      session.ls_precfns <- PrecFns { prec_solve_fn = solve;
                                      prec_setup_fn = setup }

    let prec_none = InternalPrecNone (fun session nv ->
        session.ls_precfns <- NoPrecFns)
    let prec_left ?setup solve =
      InternalPrecLeft (init_preconditioner solve setup)
    let prec_right ?setup solve =
      InternalPrecRight (init_preconditioner solve setup)
    let prec_both ?setup solve =
      InternalPrecBoth (init_preconditioner solve setup)

    let set_jac_times_vec_fn s f =
      match s.ls_callbacks with
      | SpilsCallback _ ->
          c_set_jac_times_vec_fn s true;
          s.ls_callbacks <- SpilsCallback (Some f)
      | _ -> raise Sundials.InvalidLinearSolver

    let init_spils init maxl jac_times_vec prec session nv =
      let with_prec prec_type set_prec =
        init session maxl prec_type;
        set_prec session nv;
        session.ls_callbacks <- SpilsCallback jac_times_vec;
        match jac_times_vec with
        | None -> ()
        | Some jtv -> set_jac_times_vec_fn session jtv
      in
      match prec with
      | InternalPrecNone set_prec  -> with_prec Spils.PrecNone set_prec
      | InternalPrecLeft set_prec  -> with_prec Spils.PrecLeft set_prec
      | InternalPrecRight set_prec -> with_prec Spils.PrecRight set_prec
      | InternalPrecBoth set_prec  -> with_prec Spils.PrecBoth set_prec

    let spgmr ?(maxl=0) ?jtv prec session nv =
      init_spils c_spgmr maxl jtv prec session nv

    let spbcg ?(maxl=0) ?jtv prec session nv =
      init_spils c_spbcg maxl jtv prec session nv

    let sptfqmr ?(maxl=0) ?jtv prec session nv =
      init_spils c_sptfqmr maxl jtv prec session nv

    let set_preconditioner s ?setup solve =
      match s.ls_callbacks with
      | SpilsCallback _ ->
          c_set_preconditioner s (setup <> None);
          s.ls_precfns <- PrecFns { prec_setup_fn = setup;
                                    prec_solve_fn = solve }
      | _ -> raise Sundials.InvalidLinearSolver

    let clear_jac_times_vec_fn s =
      match s.ls_callbacks with
      | SpilsCallback _ ->
          c_set_jac_times_vec_fn s false;
          s.ls_callbacks <- SpilsCallback None
      | _ -> raise Sundials.InvalidLinearSolver

    external set_prec_type
        : ('a, 'k) session -> Spils.preconditioning_type -> unit
        = "c_cvode_spils_set_prec_type"

    let set_prec_type s t =
      ls_check_spils s;
      set_prec_type s t

    external set_gs_type : ('a, 'k) session -> gramschmidt_type -> unit
        = "c_cvode_spils_set_gs_type"

    let set_gs_type s t =
      ls_check_spils s;
      set_gs_type s t

    external set_eps_lin            : ('a, 'k) session -> float -> unit
        = "c_cvode_spils_set_eps_lin"

    let set_eps_lin s epsl =
      ls_check_spils s;
      set_eps_lin s epsl

    external set_maxl               : ('a, 'k) session -> int -> unit
        = "c_cvode_spils_set_maxl"

    let set_maxl s maxl =
      ls_check_spils s;
      set_maxl s maxl

    external get_num_lin_iters      : ('a, 'k) session -> int
        = "c_cvode_spils_get_num_lin_iters"

    let get_num_lin_iters s =
      ls_check_spils s;
      get_num_lin_iters s

    external get_num_conv_fails     : ('a, 'k) session -> int
        = "c_cvode_spils_get_num_conv_fails"

    let get_num_conv_fails s =
      ls_check_spils s;
      get_num_conv_fails s

    external get_work_space         : ('a, 'k) session -> int * int
        = "c_cvode_spils_get_work_space"

    let get_work_space s =
      ls_check_spils s;
      get_work_space s

    external get_num_prec_evals     : ('a, 'k) session -> int
        = "c_cvode_spils_get_num_prec_evals"

    let get_num_prec_evals s =
      ls_check_spils s;
      get_num_prec_evals s

    external get_num_prec_solves    : ('a, 'k) session -> int
        = "c_cvode_spils_get_num_prec_solves"

    let get_num_prec_solves s =
      ls_check_spils s;
      get_num_prec_solves s

    external get_num_jtimes_evals   : ('a, 'k) session -> int
        = "c_cvode_spils_get_num_jtimes_evals"

    let get_num_jtimes_evals s =
      ls_check_spils s;
      get_num_jtimes_evals s

    external get_num_rhs_evals      : ('a, 'k) session -> int
        = "c_cvode_spils_get_num_rhs_evals"

    let get_num_rhs_evals s =
      ls_check_spils s;
      get_num_rhs_evals s

    module Banded = struct
      external c_set_preconditioner
        : ('a, 'k) session -> int -> int -> int -> unit
        = "c_cvode_spils_set_banded_preconditioner"
      (* Note: CVBandPrecInit seems to be designed only to be called on
         a fresh spils solver (i.e. right after CVSpgmr, CVSpbcg, or
         CVSptfqmr).

         In Sundials 2.5.0 and 2.6.2, calling

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

      let init_preconditioner bandrange session nv =
        c_set_preconditioner session (RealArray.length (Nvector.unwrap nv))
          bandrange.mupper bandrange.mlower;
        session.ls_precfns <- BandedPrecFns

      let prec_none = InternalPrecNone (fun session nv ->
          session.ls_precfns <- BandedPrecFns)
      let prec_left bandrange =
        InternalPrecLeft (init_preconditioner bandrange)
      let prec_right bandrange =
        InternalPrecRight (init_preconditioner bandrange)
      let prec_both bandrange =
        InternalPrecBoth (init_preconditioner bandrange)

      external get_work_space : 'k serial_session -> int * int
        = "c_cvode_bandprec_get_work_space"

      let get_work_space s =
        ls_check_spils_band s;
        get_work_space s

      external get_num_rhs_evals : 'k serial_session -> int
        = "c_cvode_bandprec_get_num_rhs_evals"

      let get_num_rhs_evals s =
        ls_check_spils_band s;
        get_num_rhs_evals s
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

    let make f s nv =
      let { linit; lsetup; lsolve } as cb = f s nv in
      c_set_alternate s (linit <> None) (lsetup <> None);
      s.ls_precfns <- NoPrecFns;
      s.ls_callbacks <- AlternateCallback cb

  end

external c_set_functional : ('a, 'k) session -> unit
  = "c_cvode_set_functional"

let set_iter_type session nv iter =
  session.ls_callbacks <- NoCallbacks;
  session.ls_precfns <- NoPrecFns;
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
  | WFtolerances of 'a error_weight_fun

let default_tolerances = SStolerances (1.0e-4, 1.0e-8)

let set_tolerances s tol =
  match tol with
  | SStolerances (rel, abs) -> (s.errw <- dummy_errw; ss_tolerances s rel abs)
  | SVtolerances (rel, abs) -> (if Sundials_config.safe then s.checkvec abs;
                                s.errw <- dummy_errw; sv_tolerances s rel abs)
  | WFtolerances ferrw -> (s.errw <- ferrw; wf_tolerances s)

external c_session_finalize : ('a, 'kind) session -> unit
    = "c_cvode_session_finalize"

let session_finalize s =
  Dls.invalidate_callback s;
  c_session_finalize s

external c_init
    : ('a, 'k) session Weak.t -> lmm -> ('a, 'k) iter -> ('a, 'k) nvector
      -> float -> (cvode_mem * c_weak_ref)
    = "c_cvode_init"

let init lmm iter tol f ?(roots=no_roots) t0 y0 =
  let (nroots, roots) = roots in
  let checkvec = Nvector.check y0 in
  if Sundials_config.safe && nroots < 0 then
    raise (Invalid_argument "number of root functions is negative");
  let weakref = Weak.create 1 in
  let cvode_mem, backref = c_init weakref lmm iter y0 t0 in
  (* cvode_mem and backref have to be immediately captured in a session and
     associated with the finalizer before we do anything else.  *)
  let session = {
          cvode        = cvode_mem;
          backref      = backref;
          nroots       = nroots;
          checkvec     = checkvec;

          exn_temp     = None;

          rhsfn        = f;
          rootsfn      = roots;
          errh         = dummy_errh;
          errw         = dummy_errw;

          ls_callbacks = NoCallbacks;
          ls_precfns   = NoPrecFns;

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

let get_num_roots { nroots } = nroots

external c_reinit
    : ('a, 'k) session -> float -> ('a, 'k) nvector -> unit
    = "c_cvode_reinit"

let reinit session ?iter ?roots t0 y0 =
  if Sundials_config.safe then session.checkvec y0;
  Dls.invalidate_callback session;
  c_reinit session t0 y0;
  (match iter with
   | None -> ()
   | Some iter_type -> set_iter_type session y0 iter_type);
  (match roots with
   | None -> ()
   | Some roots -> root_init session roots)

external get_root_info  : ('a, 'k) session -> Sundials.Roots.t -> unit
    = "c_cvode_get_root_info"

type solver_result =
  | Success             (** CV_SUCCESS *)
  | RootsFound          (** CV_ROOT_RETURN *)
  | StopTimeReached     (** CV_TSTOP_RETURN *)

external c_solve_normal : ('a, 'k) session -> float -> ('a, 'k) nvector
                              -> float * solver_result
    = "c_cvode_solve_normal"

let solve_normal s t y =
  if Sundials_config.safe then s.checkvec y;
  c_solve_normal s t y

external c_solve_one_step : ('a, 'k) session -> float -> ('a, 'k) nvector
                              -> float * solver_result
    = "c_cvode_solve_one_step"

let solve_one_step s t y =
  if Sundials_config.safe then s.checkvec y;
  c_solve_one_step s t y

external c_get_dky
    : ('a, 'k) session -> float -> int -> ('a, 'k) nvector -> unit
    = "c_cvode_get_dky"

let get_dky s y =
  if Sundials_config.safe then s.checkvec y;
  fun t k -> c_get_dky s t k y

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

let print_integrator_stats s oc =
  let stats = get_integrator_stats s
  in
    Printf.fprintf oc "num_steps = %d\n"           stats.num_steps;
    Printf.fprintf oc "num_rhs_evals = %d\n"       stats.num_rhs_evals;
    Printf.fprintf oc "num_lin_solv_setups = %d\n" stats.num_lin_solv_setups;
    Printf.fprintf oc "num_err_test_fails = %d\n"  stats.num_err_test_fails;
    Printf.fprintf oc "last_order = %d\n"          stats.last_order;
    Printf.fprintf oc "current_order = %d\n"       stats.current_order;
    Printf.fprintf oc "actual_init_step = %e\n"    stats.actual_init_step;
    Printf.fprintf oc "last_step = %e\n"           stats.last_step;
    Printf.fprintf oc "current_step = %e\n"        stats.current_step;
    Printf.fprintf oc "current_time = %e\n"        stats.current_time;

external set_error_file : ('a, 'k) session -> Sundials.Logfile.t -> unit
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
  set_root_direction' s (Sundials.RootDirs.copy (get_num_roots s) rda)

let set_all_root_directions s rd =
  set_root_direction' s (Sundials.RootDirs.make (get_num_roots s) rd)

external set_no_inactive_root_warn      : ('a, 'k) session -> unit
    = "c_cvode_set_no_inactive_root_warn"

external get_num_stab_lim_order_reds    : ('a, 'k) session -> int
    = "c_cvode_get_num_stab_lim_order_reds"

external get_tol_scale_factor           : ('a, 'k) session -> float
    = "c_cvode_get_tol_scale_factor"

external c_get_err_weights : ('a, 'k) session -> ('a, 'k) nvector -> unit
    = "c_cvode_get_err_weights"

let get_err_weights s ew =
  if Sundials_config.safe then s.checkvec ew;
  c_get_err_weights s ew

external c_get_est_local_errors : ('a, 'k) session -> ('a, 'k) nvector -> unit
    = "c_cvode_get_est_local_errors"

let get_est_local_errors s ew =
  if Sundials_config.safe then s.checkvec ew;
  c_get_est_local_errors s ew

external get_num_nonlin_solv_iters      : ('a, 'k) session -> int
    = "c_cvode_get_num_nonlin_solv_iters"

external get_num_nonlin_solv_conv_fails : ('a, 'k) session -> int
    = "c_cvode_get_num_nonlin_solv_conv_fails"

external get_nonlin_solv_stats          : ('a, 'k) session -> int * int
    = "c_cvode_get_nonlin_solv_stats"

external get_num_g_evals                : ('a, 'k) session -> int
    = "c_cvode_get_num_g_evals"


(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "c_cvode_init_module"

let _ =
  c_init_module
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
    |]
