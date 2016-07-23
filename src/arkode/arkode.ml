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

module Dls =
  struct
    include DlsTypes

    external c_dls_dense : 'k serial_session -> int -> bool -> unit
      = "c_arkode_dls_dense"

    external c_dls_lapack_dense : 'k serial_session -> int -> bool -> unit
      = "c_arkode_dls_lapack_dense"

    external c_dls_band : ('k serial_session * int) -> int -> int -> bool -> unit
      = "c_arkode_dls_band"

    external c_dls_lapack_band
      : ('k serial_session * int) -> int -> int -> bool -> unit
      = "c_arkode_dls_lapack_band"

    external set_dense_jac_fn : 'k serial_session -> unit
        = "c_arkode_dls_set_dense_jac_fn"

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
        = "c_arkode_dls_clear_dense_jac_fn"

    let clear_dense_jac_fn s =
      match s.ls_callbacks with
      | DlsDenseCallback _ ->
          invalidate_callback s;
          s.ls_callbacks <- DlsDenseCallback no_dense_callback;
          clear_dense_jac_fn s
      | _ -> raise Sundials.InvalidLinearSolver

    external set_band_jac_fn : 'k serial_session -> unit
        = "c_arkode_dls_set_band_jac_fn"

    let set_band_jac_fn s f =
      match s.ls_callbacks with
      | DlsBandCallback _ ->
          invalidate_callback s;
          s.ls_callbacks <- DlsBandCallback { bjacfn = f; bmat = None };
          set_band_jac_fn s
      | _ -> raise Sundials.InvalidLinearSolver

    external clear_band_jac_fn : 'k serial_session -> unit
        = "c_arkode_dls_clear_band_jac_fn"

    let clear_band_jac_fn s =
      match s.ls_callbacks with
      | DlsBandCallback _ ->
          invalidate_callback s;
          s.ls_callbacks <- DlsBandCallback no_band_callback;
          clear_band_jac_fn s
      | _ -> raise Sundials.InvalidLinearSolver

    external get_work_space : 'k serial_session -> int * int
        = "c_arkode_dls_get_work_space"

    let get_work_space s =
      ls_check_dls s;
      get_work_space s

    external get_num_jac_evals : 'k serial_session -> int
        = "c_arkode_dls_get_num_jac_evals"

    let get_num_jac_evals s =
      ls_check_dls s;
      get_num_jac_evals s

    external get_num_rhs_evals : 'k serial_session -> int
        = "c_arkode_dls_get_num_rhs_evals"

    let get_num_rhs_evals s =
      ls_check_dls s;
      get_num_rhs_evals s

    module Mass = struct
      include DlsTypes.MassTypes

      external c_dls_mass_dense : 'k serial_session -> int -> unit
        = "c_arkode_dls_mass_dense"

      external c_dls_mass_lapack_dense : 'k serial_session -> int -> unit
        = "c_arkode_dls_mass_lapack_dense"

      external c_dls_mass_band : ('k serial_session * int) -> int -> int -> unit
        = "c_arkode_dls_mass_band"

      external c_dls_mass_lapack_band
        : ('k serial_session * int) -> int -> int -> unit
        = "c_arkode_dls_mass_lapack_band"

      let dense f session nv =
        let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
        session.mass_callbacks
          <- DlsDenseMassCallback { massfn = f; dmat = None };
        session.mass_precfns <- NoMassPrecFns;
        c_dls_mass_dense session neqs

      let lapack_dense f session nv =
        let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
        session.mass_callbacks
          <- DlsDenseMassCallback { massfn = f; dmat = None };
        session.mass_precfns <- NoMassPrecFns;
        c_dls_mass_lapack_dense session neqs

      let band f p session nv =
        let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
        session.mass_callbacks <-
           DlsBandMassCallback { bmassfn = f; bmat = None };
        session.mass_precfns <- NoMassPrecFns;
        c_dls_mass_band (session, neqs) p.mupper p.mlower

      let lapack_band f p session nv =
        let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
        session.mass_callbacks <-
           DlsBandMassCallback { bmassfn = f; bmat = None };
        session.mass_precfns <- NoMassPrecFns;
        c_dls_mass_lapack_band (session, neqs) p.mupper p.mlower

      let invalidate_callback session =
        match session.mass_callbacks with
        | DlsDenseMassCallback ({ dmat = Some d } as cb) ->
            Dls.DenseMatrix.invalidate d;
            cb.dmat <- None
        | DlsBandMassCallback  ({ bmat = Some d } as cb) ->
            Dls.BandMatrix.invalidate d;
            cb.bmat <- None
        | SlsKluMassCallback ({ SlsTypes.smmat = Some d } as cb) ->
            Sls_impl.invalidate d;
            cb.SlsTypes.smmat <- None
        | SlsSuperlumtMassCallback ({ SlsTypes.smmat = Some d } as cb)
            -> Sls_impl.invalidate d;
               cb.SlsTypes.smmat <- None
        | _ -> ()

      let set_dense_fn s f =
        match s.mass_callbacks with
        | DlsDenseMassCallback _ ->
            invalidate_callback s;
            s.mass_callbacks <- DlsDenseMassCallback { massfn = f; dmat = None }
        | _ -> raise Sundials.InvalidLinearSolver

      let set_band_fn s f =
        match s.mass_callbacks with
        | DlsBandMassCallback _ ->
            invalidate_callback s;
            s.mass_callbacks <- DlsBandMassCallback { bmassfn = f; bmat = None }
        | _ -> raise Sundials.InvalidLinearSolver

      external get_work_space : 'k serial_session -> int * int
          = "c_arkode_dls_get_mass_work_space"

      let get_work_space s =
        mass_check_dls s;
        get_work_space s

      external get_num_evals : 'k serial_session -> int
          = "c_arkode_dls_get_num_mass_evals"

      let get_num_evals s =
        mass_check_dls s;
        get_num_evals s
    end
  end

module Sls = struct
  include SlsTypes
  
  module Klu = struct

    (* Must correspond with arkode_klu_ordering_tag *)
    type ordering =
         Amd
       | ColAmd
       | Natural

    external c_klu : 'k serial_session -> int -> int -> unit
      = "c_arkode_klu_init"

    let solver f nnz session nv =
      if not Sundials_config.klu_enabled
        then raise Sundials.NotImplementedBySundialsVersion;
      let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
      session.ls_callbacks <- SlsKluCallback { jacfn = f; smat = None };
      session.ls_precfns <- NoPrecFns;
      c_klu session neqs nnz

    external c_set_ordering : 'k serial_session -> ordering -> unit
      = "c_arkode_klu_set_ordering"

    let set_ordering session ordering =
      ls_check_klu session;
      c_set_ordering session ordering

    external c_reinit : 'k serial_session -> int -> int -> bool -> unit
      = "c_arkode_klu_reinit"

    let reinit session n nnz realloc =
      ls_check_klu session;
      c_reinit session n nnz realloc

    external c_get_num_jac_evals : 'k serial_session -> int
      = "c_arkode_klu_get_num_jac_evals"

    let get_num_jac_evals session =
      ls_check_klu session;
      c_get_num_jac_evals session

    module Mass = struct

      external c_mass_klu : 'k serial_session -> int -> int -> unit
        = "c_arkode_mass_klu_init"

      let solver f nnz session nv =
        if not Sundials_config.klu_enabled
          then raise Sundials.NotImplementedBySundialsVersion;
        let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
        session.mass_callbacks
          <- SlsKluMassCallback { massfn = f; smmat = None };
        session.mass_precfns <- NoMassPrecFns;
        c_mass_klu session neqs nnz

      external c_set_ordering : 'k serial_session -> ordering -> unit
        = "c_arkode_mass_klu_set_ordering"

      let set_ordering session ordering =
        mass_check_klu session;
        c_set_ordering session ordering

      external c_reinit : 'k serial_session -> int -> int -> bool -> unit
        = "c_arkode_mass_klu_reinit"

      let reinit session n nnz realloc =
        mass_check_klu session;
        c_reinit session n nnz realloc

      external c_get_num_evals : 'k serial_session -> int
        = "c_arkode_klu_get_num_mass_evals"

      let get_num_evals session =
        mass_check_klu session;
        c_get_num_jac_evals session
    end
  end

  module Superlumt = struct

    (* Must correspond with arkode_superlumt_ordering_tag *)
    type ordering =
         Natural
       | MinDegreeProd
       | MinDegreeSum
       | ColAmd

    external c_superlumt : 'k serial_session -> int -> int -> int -> unit
      = "c_arkode_superlumt_init"

    let solver f ~nnz ~nthreads session nv =
      if not Sundials_config.superlumt_enabled
        then raise Sundials.NotImplementedBySundialsVersion;
      let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
      session.ls_callbacks <- SlsSuperlumtCallback { jacfn = f; smat = None };
      session.ls_precfns <- NoPrecFns;
      c_superlumt session neqs nnz nthreads

    external c_set_ordering : 'k serial_session -> ordering -> unit
      = "c_arkode_superlumt_set_ordering"

    let set_ordering session ordering =
      ls_check_superlumt session;
      c_set_ordering session ordering

    external c_get_num_jac_evals : 'k serial_session -> int
      = "c_arkode_superlumt_get_num_jac_evals"

    let get_num_jac_evals session =
      ls_check_superlumt session;
      c_get_num_jac_evals session

    module Mass = struct
      external c_mass_superlumt : 'k serial_session -> int -> int -> int -> unit
        = "c_arkode_mass_superlumt_init"

      let solver f ~nnz ~nthreads session nv =
        if not Sundials_config.superlumt_enabled
          then raise Sundials.NotImplementedBySundialsVersion;
        let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
        session.mass_callbacks
          <- SlsSuperlumtMassCallback { massfn = f; smmat = None };
        session.mass_precfns <- NoMassPrecFns;
        c_mass_superlumt session neqs nnz nthreads

      external c_set_ordering : 'k serial_session -> ordering -> unit
        = "c_arkode_mass_superlumt_set_ordering"

      let set_ordering session ordering =
        mass_check_superlumt session;
        c_set_ordering session ordering

      external c_get_num_evals : 'k serial_session -> int
        = "c_arkode_superlumt_get_num_mass_evals"

      let get_num_evals session =
        mass_check_superlumt session;
        c_get_num_evals session
    end
  end
end

module Spils =
  struct
    include SpilsTypes

    external c_spgmr
      : ('a, 'k) session -> int -> Spils.preconditioning_type -> unit
      = "c_arkode_spils_spgmr"

    external c_spbcg
      : ('a, 'k) session -> int -> Spils.preconditioning_type -> unit
      = "c_arkode_spils_spbcg"

    external c_sptfqmr
      : ('a, 'k) session -> int -> Spils.preconditioning_type -> unit
      = "c_arkode_spils_sptfqmr"

    external c_spfgmr
      : ('a, 'k) session -> int -> Spils.preconditioning_type -> unit
      = "c_arkode_spils_spfgmr"

    external c_pcg
      : ('a, 'k) session -> int -> Spils.preconditioning_type -> unit
      = "c_arkode_spils_pcg"

    external c_set_preconditioner
      : ('a, 'k) session -> bool -> unit
      = "c_arkode_spils_set_preconditioner"

    external c_set_jac_times_vec_fn
      : ('a, 'k) session -> bool -> unit
      = "c_arkode_spils_set_jac_times_vec_fn"

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
      | SpilsCallback cbs ->
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
      | InternalPrecNone set_prec  -> with_prec Spils.PrecNone set_prec;
      | InternalPrecLeft set_prec  -> with_prec Spils.PrecLeft set_prec
      | InternalPrecRight set_prec -> with_prec Spils.PrecRight set_prec
      | InternalPrecBoth set_prec  -> with_prec Spils.PrecBoth set_prec

    let spgmr ?(maxl=0) ?jac_times_vec prec session nv =
      init_spils c_spgmr maxl jac_times_vec prec session nv

    let spbcg ?(maxl=0) ?jac_times_vec prec session nv =
      init_spils c_spbcg maxl jac_times_vec prec session nv

    let sptfqmr ?(maxl=0) ?jac_times_vec prec session nv =
      init_spils c_sptfqmr maxl jac_times_vec prec session nv

    let spfgmr ?(maxl=0) ?jac_times_vec prec session nv =
      init_spils c_spfgmr maxl jac_times_vec prec session nv

    let pcg ?(maxl=0) ?jac_times_vec prec session nv =
      init_spils c_pcg maxl jac_times_vec prec session nv

    let set_preconditioner s ?setup solve =
      match s.ls_callbacks with
      | SpilsCallback _ ->
          c_set_preconditioner s (setup <> None);
          s.ls_precfns <- PrecFns { prec_setup_fn = setup;
                                    prec_solve_fn = solve }
      | _ -> raise Sundials.InvalidLinearSolver

    let clear_jac_times_vec_fn s =
      match s.ls_callbacks with
      | SpilsCallback cbs ->
          c_set_jac_times_vec_fn s false;
          s.ls_callbacks <- SpilsCallback None
      | _ -> raise Sundials.InvalidLinearSolver

    external set_prec_type
        : ('a, 'k) session -> Spils.preconditioning_type -> unit
        = "c_arkode_spils_set_prec_type"

    let set_prec_type s t =
      ls_check_spils s;
      set_prec_type s t

    external set_gs_type : ('a, 'k) session -> gramschmidt_type -> unit
        = "c_arkode_spils_set_gs_type"

    let set_gs_type s t =
      ls_check_spils s;
      set_gs_type s t

    external set_eps_lin            : ('a, 'k) session -> float -> unit
        = "c_arkode_spils_set_eps_lin"

    let set_eps_lin s epsl =
      ls_check_spils s;
      set_eps_lin s epsl

    external set_maxl               : ('a, 'k) session -> int -> unit
        = "c_arkode_spils_set_maxl"

    let set_maxl s maxl =
      ls_check_spils s;
      set_maxl s maxl

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

    module Banded = struct
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
        = "c_arkode_bandprec_get_work_space"

      let get_work_space s =
        ls_check_spils_band s;
        get_work_space s

      external get_num_rhs_evals : 'k serial_session -> int
        = "c_arkode_bandprec_get_num_rhs_evals"

      let get_num_rhs_evals s =
        ls_check_spils_band s;
        get_num_rhs_evals s
    end

    module Mass = struct
      include SpilsTypes.MassTypes

      external c_spgmr
        : ('a, 'k) session -> int -> Spils.preconditioning_type -> unit
        = "c_arkode_spils_mass_spgmr"

      external c_spbcg
        : ('a, 'k) session -> int -> Spils.preconditioning_type -> unit
        = "c_arkode_spils_mass_spbcg"

      external c_sptfqmr
        : ('a, 'k) session -> int -> Spils.preconditioning_type -> unit
        = "c_arkode_spils_mass_sptfqmr"

      external c_spfgmr
        : ('a, 'k) session -> int -> Spils.preconditioning_type -> unit
        = "c_arkode_spils_mass_spfgmr"

      external c_pcg
        : ('a, 'k) session -> int -> Spils.preconditioning_type -> unit
        = "c_arkode_spils_mass_pcg"

      external c_set_preconditioner
        : ('a, 'k) session -> bool -> unit
        = "c_arkode_spils_set_mass_preconditioner"

      let init_preconditioner solve setup session nv =
        c_set_preconditioner session (setup <> None);
        session.mass_precfns <- MassPrecFns { prec_solve_fn = solve;
                                              prec_setup_fn = setup }

      let prec_none = InternalPrecNone (fun session nv ->
          session.mass_precfns <- NoMassPrecFns)
      let prec_left ?setup solve =
        InternalPrecLeft (init_preconditioner solve setup)
      let prec_right ?setup solve =
        InternalPrecRight (init_preconditioner solve setup)
      let prec_both ?setup solve =
        InternalPrecBoth (init_preconditioner solve setup)

      let set_times_vec_fn s f =
        match s.mass_callbacks with
        | SpilsMassCallback _ ->
            s.mass_callbacks <- SpilsMassCallback f
        | _ -> raise Sundials.InvalidLinearSolver

      let init_spils init maxl times prec session nv =
        let with_prec prec_type set_prec =
          init session maxl prec_type;
          set_prec session nv;
          session.mass_callbacks <- SpilsMassCallback times;
          set_times_vec_fn session times
        in
        match prec with
        | InternalPrecNone set_prec  -> with_prec Spils.PrecNone  set_prec
        | InternalPrecLeft set_prec  -> with_prec Spils.PrecLeft  set_prec
        | InternalPrecRight set_prec -> with_prec Spils.PrecRight set_prec
        | InternalPrecBoth set_prec  -> with_prec Spils.PrecBoth  set_prec

      let spgmr ?(maxl=0) times prec session nv =
        init_spils c_spgmr maxl times prec session nv

      let spbcg ?(maxl=0) times prec session nv =
        init_spils c_spbcg maxl times prec session nv

      let sptfqmr ?(maxl=0) times prec session nv =
        init_spils c_sptfqmr maxl times prec session nv

      let spfgmr ?(maxl=0) times prec session nv =
        init_spils c_spfgmr maxl times prec session nv

      let pcg ?(maxl=0) times prec session nv =
        init_spils c_pcg maxl times prec session nv

      let set_preconditioner s ?setup solve =
        match s.mass_callbacks with
        | SpilsMassCallback _ ->
            c_set_preconditioner s (setup <> None);
            s.mass_precfns <- MassPrecFns { prec_setup_fn = setup;
                                            prec_solve_fn = solve }
        | _ -> raise Sundials.InvalidLinearSolver

      external set_prec_type
          : ('a, 'k) session -> Spils.preconditioning_type -> unit
          = "c_arkode_spils_set_mass_prec_type"

      let set_prec_type s t =
        mass_check_spils s;
        set_prec_type s t

      external set_gs_type : ('a, 'k) session -> gramschmidt_type -> unit
          = "c_arkode_spils_set_mass_gs_type"

      let set_gs_type s t =
        mass_check_spils s;
        set_gs_type s t

      external set_eps_lin            : ('a, 'k) session -> float -> unit
          = "c_arkode_spils_set_mass_eps_lin"

      let set_eps_lin s epsl =
        mass_check_spils s;
        set_eps_lin s epsl

      external set_maxl               : ('a, 'k) session -> int -> unit
          = "c_arkode_spils_set_mass_maxl"

      let set_maxl s maxl =
        mass_check_spils s;
        set_maxl s maxl

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
    end
  end

module Alternate =
  struct
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

    module Mass = struct
      include AlternateTypes.MassTypes

      external c_set_alternate
        : ('data, 'kind) session -> bool -> bool -> unit
        = "c_arkode_set_mass_alternate"

      let make f s nv =
        let { minit; msetup; msolve } as cb = f s nv in
        c_set_alternate s (minit <> None) (msetup <> None);
        s.mass_precfns <- NoMassPrecFns;
        s.mass_callbacks <- AlternateMassCallback cb
    end
  end

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
  -> (arkode_mem * c_weak_ref * arkode_file)
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
  let arkode_mem, backref, no_file
        = c_init weakref (fi <> None) (fe <> None) y0 t0 in
  (* arkode_mem and backref have to be immediately captured in a session and
     associated with the finalizer before we do anything else.  *)
  let session = {
          arkode       = arkode_mem;
          backref      = backref;
          nroots       = nroots;
          err_file     = no_file;
          diag_file    = no_file;
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

          linsolver      = None;
          ls_callbacks   = NoCallbacks;
          ls_precfns     = NoPrecFns;
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

external set_diagnostics : ('a, 'k) session -> string -> bool -> unit
    = "c_arkode_set_diagnostics"

external set_error_file : ('a, 'k) session -> string -> bool -> unit
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
    stage_times : RealArray.t;
    coefficients : RealArray.t;
    bembed : RealArray.t;
  }

external c_set_ark_tables
  : ('d, 'k) session -> rk_method -> RealArray.t -> RealArray.t -> unit
    = "c_arkode_set_ark_tables"

external c_set_erk_table
  : ('d, 'k) session -> rk_method -> RealArray.t -> unit
    = "c_arkode_set_erk_table"

external c_set_irk_table
  : ('d, 'k) session -> rk_method -> RealArray.t -> unit
    = "c_arkode_set_irk_table"

let set_ark_tables s rkm ai ae =
  (if s.irhsfn == dummy_irhsfn || s.erhsfn == dummy_erhsfn then raise IllInput);
  c_set_ark_tables s rkm ai ae

let set_erk_table s rkm ae =
  (if s.erhsfn == dummy_erhsfn then raise IllInput);
  c_set_erk_table s rkm ae

let set_irk_table s rkm ai =
  (if s.irhsfn == dummy_irhsfn then raise IllInput);
  c_set_irk_table s rkm ai

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
  match v with
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

let int_of_irk_table v =
  match v with
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
  | ARK_8_4_5_Implicit -> 22

let ints_of_ark_table v =
  match v with
  | ARK_4_2_3 -> (15, 2)
  | ARK_6_3_4 -> (20, 4)
  | ARK_8_4_5 -> (22, 9)

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

external c_set_root_direction   : ('a, 'k) session -> Sundials.RootDirs.t -> unit
    = "c_arkode_set_root_direction"

let set_root_direction s rda =
  c_set_root_direction s (Sundials.RootDirs.copy (get_num_roots s) rda)

let set_all_root_directions s rd =
  c_set_root_direction s (Sundials.RootDirs.make (get_num_roots s) rd)

external set_no_inactive_root_warn      : ('a, 'k) session -> unit
    = "c_arkode_set_no_inactive_root_warn"

external get_current_butcher_tables
  : ('d, 'k) session -> RealArray.t * RealArray.t * rk_method
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
      BadK;
      BadT;
    |]
