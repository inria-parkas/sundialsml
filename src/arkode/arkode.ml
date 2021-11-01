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
open Sundials

(*
 * NB: The order of variant constructors and record fields is important!
 *     If these types are changed or augmented, the corresponding declarations
 *     in arkode_serial.h (and code in arkode_serial.c) must also be updated.
 *)

(* Solver exceptions *)
exception IllInput
exception TooClose
exception TooMuchWork
exception TooMuchAccuracy
exception InnerStepFail of exn option
exception ErrFailure
exception ConvergenceFailure
exception LinearInitFailure
exception LinearSetupFailure of exn option
exception LinearSolveFailure of exn option
exception NonlinearInitFailure
exception NonlinearSetupFailure
exception NonlinearSetupRecoverable
exception NonlinearOperationError
exception MassInitFailure
exception MassSetupFailure of exn option
exception MassSolveFailure of exn option
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

exception VectorOpErr

type arkstep = Arkode_impl.arkstep
type erkstep = Arkode_impl.erkstep
type mristep = Arkode_impl.mristep

module Common = struct (* {{{ *)

  include Arkode_impl.Global

  let no_roots = (0, Arkode_impl.dummy_rootsfn)

  (* Synchronized with arkode_step_stats_index in arkode_ml.h *)
  type step_stats = {
      num_steps           : int;
      actual_init_step    : float;
      last_step           : float;
      current_step        : float;
      current_time        : float
    }

  let print_step_stats oc stats =
    Printf.fprintf oc "num_steps = %d\n"           stats.num_steps;
    Printf.fprintf oc "actual_init_step = %e\n"    stats.actual_init_step;
    Printf.fprintf oc "last_step = %e\n"           stats.last_step;
    Printf.fprintf oc "current_step = %e\n"        stats.current_step;
    Printf.fprintf oc "current_time = %e\n"        stats.current_time

  (* Synchronized with arkode_solver_result_tag in arkode_ml.h *)
  type solver_result =
    | Success             (** ARK_SUCCESS *)
    | RootsFound          (** ARK_ROOT_RETURN *)
    | StopTimeReached     (** ARK_TSTOP_RETURN *)

  (* must correspond with the static table ark_interpolant_types in arkode_ml.c *)
  type interpolant_type =
      Hermite
    | Lagrange

  type linearity =
    | Linear of bool
    | Nonlinear

  type ('a, 'k) tolerance =
    | SStolerances of float * float
    | SVtolerances of float * ('a, 'k) Nvector.nvector
    | WFtolerances of 'a error_weight_fun

  let default_tolerances = SStolerances (1.0e-4, 1.0e-9)

  type predictor_method =
    | TrivialPredictor
    | MaximumOrderPredictor
    | VariableOrderPredictor
    | CutoffOrderPredictor
    | BootstrapPredictor
    | MinimumCorrectionPredictor

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

  (* must correspond to arkode_nonlin_system_data_index in arkode_ml.h *)
  type 'd nonlin_system_data = {
    tcur  : float;
    zpred : 'd;
    zi    : 'd;
    fi    : 'd;
    gamma : float;
    sdata : 'd;
  }

  type 'd triple = 'd * 'd * 'd

  type ('t, 'd) jacobian_arg = ('t, 'd) Arkode_impl.jacobian_arg =
    {
      jac_t   : float;
      jac_y   : 'd;
      jac_fy  : 'd;
      jac_tmp : 't;
    }

end (* }}} *)

module Dls = struct
  include Arkode_impl.DirectTypes
end

module Spils = struct
  include Arkode_impl.SpilsTypes
  module LSI = Sundials_LinearSolver_impl

  module Banded = struct (* {{{ *)
    (* These fields are accessed from arkode_ml.c *)
    type bandrange = { mupper : int; mlower : int; }

    external c_set_preconditioner
      : ('a, 'k, 's) Arkode_impl.session -> int -> int -> int -> unit
      = "sunml_arkode_set_banded_preconditioner"
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
      LSI.Iterative.(PrecNone,
                     fun session nv ->
                            Arkode_impl.(session.ls_precfns <- BandedPrecFns))
    let prec_left bandrange =
      LSI.Iterative.(PrecLeft,  init_preconditioner bandrange)
    let prec_right bandrange =
      LSI.Iterative.(PrecRight, init_preconditioner bandrange)
    let prec_both bandrange =
      LSI.Iterative.(PrecBoth,  init_preconditioner bandrange)

    type ('k, 's) serial_session =
      (Nvector_serial.data, 'k, 's) Arkode_impl.session
      constraint 'k = [>Nvector_serial.kind]
      constraint 's = [<arkstep|mristep]

    external get_work_space : ('k, 's) Arkode_impl.serial_session -> int * int
      = "sunml_arkode_bandprec_get_work_space"

    let get_work_space s =
      Arkode_impl.ls_check_spils_band s;
      get_work_space s

    external get_num_rhs_evals : ('k, 's) Arkode_impl.serial_session -> int
      = "sunml_arkode_bandprec_get_num_rhs_evals"

    let get_num_rhs_evals s =
      Arkode_impl.ls_check_spils_band s;
      get_num_rhs_evals s
  end (* }}} *)

end

module ButcherTable = struct (* {{{ *)

  (* Synchronized with arkode_butcher_table_index in arkode_ml.h *)
  type t = {
      method_order : int;
      embedding_order : int;
      stages : int;
      stage_values : Sundials.RealArray2.t;
      stage_times : RealArray.t;
      coefficients : RealArray.t;
      bembed : RealArray.t option;
    }

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
    | Fehlberg_13_7_8     (* >= 2.7.0 *)
    | Knoth_Wolke_3_3     (* >= 4.0.0 *)
    | ARK_7_3_4_Explicit  (* >= 5.0.0 *)
    | ARK_8_4_5b_Explicit (* >= 5.0.0 *)

  type dirk_table =
    | SDIRK_2_1_2
    | Billington_3_2_3
    | TRBDF2_3_2_3
    | Kvaerno_4_2_3
    | ARK324L2SA_DIRK_4_2_3
    | Cash_5_2_4
    | Cash_5_3_4
    | SDIRK_5_3_4
    | Kvaerno_5_3_4
    | ARK_6_3_4_Implicit
    | Kvaerno_7_4_5
    | ARK_8_4_5_Implicit
    | ARK_7_3_4_Implicit
    | ARK_8_4_5b_Implicit

  type ark_table =
    | ARK_4_2_3
    | ARK_6_3_4
    | ARK_8_4_5

  let int_of_erk_table v =
    match Config.sundials_version with
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
       | Fehlberg_13_7_8
       | _                     -> raise Config.NotImplementedBySundialsVersion)
    | 2,7,_ | 3,_,_ ->
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
       | Fehlberg_13_7_8       -> 11
       | _                     -> raise Config.NotImplementedBySundialsVersion)
    | 4,_,_ ->
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
       | Fehlberg_13_7_8       -> 11
       | Knoth_Wolke_3_3       -> 12
       | _                     -> raise Config.NotImplementedBySundialsVersion)
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
       | Fehlberg_13_7_8       -> 11
       | Knoth_Wolke_3_3       -> 12
       | ARK_7_3_4_Explicit    -> 13
       | ARK_8_4_5b_Explicit   -> 14)

  let int_of_dirk_table v =
    match Config.sundials_version with
    | 2,5,_ | 2,6,_ ->
      (match v with
       | SDIRK_2_1_2        -> 11
       | Billington_3_2_3   -> 12
       | TRBDF2_3_2_3       -> 13
       | Kvaerno_4_2_3      -> 14
       | ARK324L2SA_DIRK_4_2_3 -> 15
       | Cash_5_2_4         -> 16
       | Cash_5_3_4         -> 17
       | SDIRK_5_3_4        -> 18
       | Kvaerno_5_3_4      -> 19
       | ARK_6_3_4_Implicit -> 20
       | Kvaerno_7_4_5      -> 21
       | ARK_8_4_5_Implicit -> 22
       | _                  -> raise Config.NotImplementedBySundialsVersion)
    | 2,7,_ | 3,_,_ ->
      (match v with
       | SDIRK_2_1_2        -> 12
       | Billington_3_2_3   -> 13
       | TRBDF2_3_2_3       -> 14
       | Kvaerno_4_2_3      -> 15
       | ARK324L2SA_DIRK_4_2_3 -> 16
       | Cash_5_2_4         -> 17
       | Cash_5_3_4         -> 18
       | SDIRK_5_3_4        -> 19
       | Kvaerno_5_3_4      -> 20
       | ARK_6_3_4_Implicit -> 21
       | Kvaerno_7_4_5      -> 22
       | ARK_8_4_5_Implicit -> 23
       | _                  -> raise Config.NotImplementedBySundialsVersion)
    | 4,_,_ ->
      (match v with
       | SDIRK_2_1_2        -> 100
       | Billington_3_2_3   -> 101
       | TRBDF2_3_2_3       -> 102
       | Kvaerno_4_2_3      -> 103
       | ARK324L2SA_DIRK_4_2_3 -> 104
       | Cash_5_2_4         -> 105
       | Cash_5_3_4         -> 106
       | SDIRK_5_3_4        -> 107
       | Kvaerno_5_3_4      -> 108
       | ARK_6_3_4_Implicit -> 109
       | Kvaerno_7_4_5      -> 110
       | ARK_8_4_5_Implicit -> 111
       | _                  -> raise Config.NotImplementedBySundialsVersion)
    | _ ->
      (match v with
       | SDIRK_2_1_2        -> 100
       | Billington_3_2_3   -> 101
       | TRBDF2_3_2_3       -> 102
       | Kvaerno_4_2_3      -> 103
       | ARK324L2SA_DIRK_4_2_3 -> 104
       | Cash_5_2_4         -> 105
       | Cash_5_3_4         -> 106
       | SDIRK_5_3_4        -> 107
       | Kvaerno_5_3_4      -> 108
       | ARK_6_3_4_Implicit -> 109
       | Kvaerno_7_4_5      -> 110
       | ARK_8_4_5_Implicit -> 111
       | ARK_7_3_4_Implicit -> 112
       | ARK_8_4_5b_Implicit -> 113)

  let ints_of_ark_table v =
    match v with
    | ARK_4_2_3 -> (int_of_dirk_table ARK324L2SA_DIRK_4_2_3,
                    int_of_erk_table ARK_4_2_3_Explicit)
    | ARK_6_3_4 -> (int_of_dirk_table ARK_6_3_4_Implicit,
                    int_of_erk_table ARK_6_3_4_Explicit)
    | ARK_8_4_5 -> (int_of_dirk_table ARK_8_4_5_Implicit,
                    int_of_erk_table ARK_8_4_5_Explicit)

  external c_load_erk : int -> t
    = "sunml_arkode_butcher_table_load_erk"

  external c_load_dirk : int -> t
    = "sunml_arkode_butcher_table_load_dirk"

  let load_erk v = c_load_erk (int_of_erk_table v)
  let load_dirk v = c_load_dirk (int_of_dirk_table v)

  external write : t -> Logfile.t -> unit
    = "sunml_arkode_butcher_table_write"

  exception ButcherTableCheckFailed

  external c_check_order
    : Logfile.t option -> t -> int * int option * bool * bool
    = "sunml_arkode_butcher_table_check_order"

  let check_order ?outfile bt =
    let q, op, err, warn = c_check_order outfile bt in
    if err then raise ButcherTableCheckFailed;
    q, op, warn

  external c_check_ark_order
    : Logfile.t option -> t -> t -> int * int option * bool
    = "sunml_arkode_butcher_table_check_ark_order"

  let check_ark_order ?outfile = c_check_ark_order outfile

end (* }}} *)

module ARKStep = struct (* {{{ *)
  include Arkode_impl
  include Common

  type ('d, 'k) session = ('d, 'k, arkstep) Arkode_impl.session

  type 'k serial_session = (Nvector_serial.data, 'k) session
                           constraint 'k = [>Nvector_serial.kind]

  type ('d, 'k) linear_solver = ('d, 'k, arkstep) lin_solver

  type 'k serial_linear_solver =
    (Nvector_serial.data, 'k) linear_solver
    constraint 'k = [>Nvector_serial.kind]

  type ('d, 'k) implicit_problem = {
      irhsfn    : 'd rhsfn;
      linearity : linearity;
      nlsolver  : ('d, 'k, ('d, 'k) session, [`Nvec]) Sundials_NonlinearSolver.t option;
      lsolver   : ('d, 'k) linear_solver option;
    }

  let implicit_problem ?nlsolver ?lsolver ?(linearity=Nonlinear) fi =
    { irhsfn = fi; linearity; nlsolver; lsolver }

  type ('d, 'k) problem =
    | Implicit of ('d, 'k) implicit_problem
    | Explicit of 'd rhsfn
    | ImEx of 'd rhsfn * ('d, 'k) implicit_problem

  let implicit ?nlsolver ?lsolver ?linearity fi =
    Implicit (implicit_problem ?nlsolver ?lsolver ?linearity fi)

  let explicit f = Explicit f

  let imex ?nlsolver ?lsolver ?linearity ~fi fe =
    ImEx (fe, implicit_problem ?nlsolver ?lsolver ?linearity fi)

  external c_root_init : ('a, 'k) session -> int -> unit
      = "sunml_arkode_ark_root_init"

  let root_init session (nroots, rootsfn) =
    c_root_init session nroots;
    session.rootsfn <- rootsfn

  module Dls = struct (* {{{ *)
    include Dls
    include LinearSolver.Direct

    (* Sundials < 3.0.0 *)
    external c_dls_dense : 'k serial_session -> int -> bool -> unit
      = "sunml_arkode_dls_dense"

    (* Sundials < 3.0.0 *)
    external c_dls_lapack_dense
      : 'k serial_session -> int -> bool -> unit
      = "sunml_arkode_dls_lapack_dense"

    (* Sundials < 3.0.0 *)
    external c_dls_band
      : 'k serial_session -> int -> int -> int -> bool -> unit
      = "sunml_arkode_dls_band"

    (* Sundials < 3.0.0 *)
    external c_dls_lapack_band
      : 'k serial_session -> int -> int -> int -> bool -> unit
      = "sunml_arkode_dls_lapack_band"

    (* Sundials < 3.0.0 *)
    external c_klu
      : 'k serial_session -> 's Matrix.Sparse.sformat
        -> int -> int -> unit
      = "sunml_arkode_klu_init"

    (* Sundials < 3.0.0 *)
    external c_klu_set_ordering
      : 'k serial_session -> LinearSolver.Direct.Klu.ordering -> unit
      = "sunml_arkode_klu_set_ordering"

    (* Sundials < 3.0.0 *)
    external c_klu_reinit
      : 'k serial_session -> int -> int -> unit
      = "sunml_arkode_klu_reinit"

    (* Sundials < 3.0.0 *)
    external c_superlumt
      : 'k serial_session -> int -> int -> int -> unit
      = "sunml_arkode_superlumt_init"

    (* Sundials < 3.0.0 *)
    external c_superlumt_set_ordering
      : 'k serial_session -> LinearSolver.Direct.Superlumt.ordering
        -> unit
      = "sunml_arkode_superlumt_set_ordering"

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
          (solver_data : (s, 'nd, 'nk, tag) LSI.solver_data)
          (mat : ('k, s, 'nd, 'nk) Matrix.t) session =
      match solver_data with
      | LSI.Dense ->
          let m, n = Matrix.(Dense.size (unwrap mat)) in
          if m <> n then raise LinearSolver.MatrixNotSquare;
          c_dls_dense session m hasjac
      | LSI.LapackDense ->
          let m, n = Matrix.(Dense.size (unwrap mat)) in
          if m <> n then raise LinearSolver.MatrixNotSquare;
          c_dls_lapack_dense session m hasjac

      | LSI.Band ->
          let open Matrix.Band in
          let { n; mu; ml } = dims (Matrix.unwrap mat) in
          c_dls_band session n mu ml hasjac
      | LSI.LapackBand ->
          let open Matrix.Band in
          let { n; mu; ml } = dims (Matrix.unwrap mat) in
          c_dls_lapack_band session n mu ml hasjac

      | LSI.Klu sinfo ->
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

      | LSI.Superlumt sinfo ->
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

      | _ -> assert false

    let check_dqjac (type k m nd nk) jac (mat : (k,m,nd,nk) Matrix.t) =
      let open Matrix in
      match get_id mat with
      | Dense | Band -> ()
      | _ -> if jac = None then invalid_arg "A Jacobian function is required"

    let set_ls_callbacks (type m) (type tag)
          ?jac ?(linsys : m linsys_fn option)
          (solver_data : (m, 'nd, 'nk, tag) LSI.solver_data)
          (mat : ('mk, m, 'nd, 'nk) Matrix.t) session =
      let cb = { jacfn = (match jac with None -> no_callback | Some f -> f);
                 jmat  = (None : m option) } in
      let ls = match linsys with None -> DirectTypes.no_linsysfn | Some f -> f in
      begin match solver_data with
      | LSI.Dense ->
          session.ls_callbacks <- DlsDenseCallback (cb, ls)
      | LSI.LapackDense ->
          session.ls_callbacks <- DlsDenseCallback (cb, ls)
      | LSI.Band ->
          session.ls_callbacks <- DlsBandCallback (cb, ls)
      | LSI.LapackBand ->
          session.ls_callbacks <- DlsBandCallback (cb, ls)
      | LSI.Klu _ ->
          if jac = None then invalid_arg "Klu requires Jacobian function";
          session.ls_callbacks <- SlsKluCallback (cb, ls)
      | LSI.Superlumt _ ->
          if jac = None then invalid_arg "Superlumt requires Jacobian function";
          session.ls_callbacks <- SlsSuperlumtCallback (cb, ls)
      | LSI.Custom _ ->
          check_dqjac jac mat;
          session.ls_callbacks <- DirectCustomCallback (cb, ls)
      | _ -> assert false
      end;
      session.ls_precfns <- NoPrecFns

    (* 3.0.0 <= Sundials < 4.0.0 *)
    external c_dls_set_linear_solver
      : 'k serial_session
        -> ('m, Nvector_serial.data, 'k) LSI.cptr
        -> ('mk, 'm, Nvector_serial.data, 'k) Matrix.t
        -> bool
        -> unit
      = "sunml_arkode_dls_set_linear_solver"

    (* 4.0.0 <= Sundials *)
    external c_set_linear_solver
      : ('d, 'k) session
        -> ('m, 'd, 'k) LSI.cptr
        -> ('mk, 'm, 'd, 'k) Matrix.t option
        -> bool
        -> bool
        -> unit
      = "sunml_arkode_ark_set_linear_solver"

    let assert_matrix = function
      | Some m -> m
      | None -> failwith "a direct linear solver is required"

    let solver ?jac ?linsys ls session nv =
      let LSI.LS ({ rawptr; solver; matrix } as hls) = ls in
      let matrix = assert_matrix matrix in
      if Sundials_impl.Version.lt500 && linsys <> None
        then raise Config.NotImplementedBySundialsVersion;
      set_ls_callbacks ?jac ?linsys solver matrix session;
      if Sundials_impl.Version.in_compat_mode2
         then make_compat (jac <> None) solver matrix session
      else if Sundials_impl.Version.in_compat_mode2_3
           then c_dls_set_linear_solver session rawptr matrix (jac <> None)
      else c_set_linear_solver session rawptr (Some matrix) (jac <> None)
                                                            (linsys <> None);
      LSI.attach ls;
      session.ls_solver <- LSI.HLS hls

    (* Sundials < 3.0.0 *)
    let invalidate_callback session =
      if Sundials_impl.Version.in_compat_mode2 then
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
        = "sunml_arkode_ark_get_lin_work_space"

    let get_work_space s =
      if Sundials_impl.Version.in_compat_mode2_3 then ls_check_direct s;
      get_work_space s

    external c_get_num_jac_evals : 'k serial_session -> int
      = "sunml_arkode_ark_get_num_jac_evals"

    (* Sundials < 3.0.0 *)
    external c_klu_get_num_jac_evals : 'k serial_session -> int
      = "sunml_arkode_klu_get_num_jac_evals"

    (* Sundials < 3.0.0 *)
    external c_superlumt_get_num_jac_evals : 'k serial_session -> int
      = "sunml_arkode_superlumt_get_num_jac_evals"

    let compat_get_num_jac_evals s =
      match s.ls_callbacks with
      | SlsKluCallback _ -> c_klu_get_num_jac_evals s
      | SlsSuperlumtCallback _ -> c_superlumt_get_num_jac_evals s
      | _ -> c_get_num_jac_evals s

    let get_num_jac_evals s =
      if Sundials_impl.Version.in_compat_mode2_3 then ls_check_direct s;
      if Sundials_impl.Version.in_compat_mode2 then compat_get_num_jac_evals s else
      c_get_num_jac_evals s

    external get_num_lin_rhs_evals : 'k serial_session -> int
        = "sunml_arkode_dls_get_num_lin_rhs_evals"

    let get_num_lin_rhs_evals s =
      if Sundials_impl.Version.in_compat_mode2_3 then ls_check_direct s;
      get_num_lin_rhs_evals s

  end (* }}} *)

  module Spils = struct (* {{{ *)
    include Spils
    include LinearSolver.Iterative

    (* Sundials < 3.0.0 *)
    external c_spgmr
      : ('a, 'k) session
        -> int -> LSI.Iterative.preconditioning_type -> unit
      = "sunml_arkode_spils_spgmr"

    (* Sundials < 3.0.0 *)
    external c_spbcgs
      : ('a, 'k) session
        -> int -> LSI.Iterative.preconditioning_type -> unit
      = "sunml_arkode_spils_spbcgs"

    (* Sundials < 3.0.0 *)
    external c_sptfqmr
      : ('a, 'k) session
        -> int -> LSI.Iterative.preconditioning_type -> unit
      = "sunml_arkode_spils_sptfqmr"

    (* Sundials < 3.0.0 *)
    external c_spfgmr
      : ('a, 'k) session
        -> int -> LSI.Iterative.preconditioning_type -> unit
      = "sunml_arkode_spils_spfgmr"

    (* Sundials < 3.0.0 *)
    external c_pcg
      : ('a, 'k) session
        -> int -> LSI.Iterative.preconditioning_type -> unit
      = "sunml_arkode_spils_pcg"

    (* Sundials < 3.0.0 *)
    external c_set_gs_type
      : ('a, 'k) session -> LSI.Iterative.gramschmidt_type -> unit
      = "sunml_arkode_spils_set_gs_type"

    (* Sundials < 3.0.0 *)
    external c_set_maxl
      : ('a, 'k) session -> int -> unit
      = "sunml_arkode_spils_set_maxl"

    (* Sundials < 3.0.0 *)
    external c_set_prec_type
      : ('a, 'k) session -> LSI.Iterative.preconditioning_type -> unit
      = "sunml_arkode_spils_set_prec_type"

    let old_set_maxl s maxl =
      if Sundials_impl.Version.in_compat_mode2_3 then ls_check_spils s;
      c_set_maxl s maxl

    let old_set_prec_type s t =
      if Sundials_impl.Version.in_compat_mode2_3 then ls_check_spils s;
      c_set_prec_type s t

    let old_set_gs_type s t =
      if Sundials_impl.Version.in_compat_mode2_3 then ls_check_spils s;
      c_set_gs_type s t

    external c_set_jac_times : ('a, 'k) session -> bool -> bool -> unit
      = "sunml_arkode_ark_set_jac_times"

    external c_set_jac_times_rhsfn : ('a, 'k) session -> bool -> unit
      = "sunml_arkode_ark_set_jac_times_rhsfn"

    external c_set_preconditioner
      : ('a, 'k) session -> bool -> unit
      = "sunml_arkode_ark_set_preconditioner"

    (* Sundials < 4.0.0 *)
    external c_spils_set_linear_solver
      : ('a, 'k) session -> ('m, 'a, 'k) LSI.cptr -> unit
      = "sunml_arkode_spils_set_linear_solver"

    (* 4.0.0 <= Sundials *)
    external c_set_linear_solver
      : ('d, 'k) session
        -> ('m, 'd, 'k) LSI.cptr
        -> ('mk, 'm, 'd, 'k) Matrix.t option
        -> bool
        -> bool
        -> unit
      = "sunml_arkode_ark_set_linear_solver"

    let init_preconditioner solve setup session nv =
      c_set_preconditioner session (setup <> None);
      session.ls_precfns <- PrecFns { prec_solve_fn = solve;
                                      prec_setup_fn = setup }

    let prec_none = LSI.Iterative.(PrecNone,
                      fun session nv -> session.ls_precfns <- NoPrecFns)

    let prec_left ?setup solve  = LSI.Iterative.(PrecLeft,
                                              init_preconditioner solve setup)

    let prec_right ?setup solve = LSI.Iterative.(PrecRight,
                                              init_preconditioner solve setup)

    let prec_both ?setup solve  = LSI.Iterative.(PrecBoth,
                                              init_preconditioner solve setup)

    (* Sundials < 3.0.0 *)
    let make_compat (type tag)
          (LSI.Iterative.({ maxl; gs_type }) as compat)
          prec_type
          (solver_data : ('s, 'nd, 'nk, tag) LSI.solver_data) session =
      match solver_data with
      | LSI.Spgmr ->
          c_spgmr session maxl prec_type;
          (match gs_type with None -> () | Some t -> c_set_gs_type session t);
          compat.set_gs_type <- old_set_gs_type session;
          compat.set_prec_type <- old_set_prec_type session
      | LSI.Spfgmr ->
          c_spfgmr session maxl prec_type;
          (match gs_type with None -> () | Some t -> c_set_gs_type session t);
          compat.set_gs_type <- old_set_gs_type session;
          compat.set_prec_type <- old_set_prec_type session
      | LSI.Spbcgs ->
          c_spbcgs session maxl prec_type;
          compat.set_maxl <- old_set_maxl session;
          compat.set_prec_type <- old_set_prec_type session
      | LSI.Sptfqmr ->
          c_sptfqmr session maxl prec_type;
          compat.set_maxl <- old_set_maxl session;
          compat.set_prec_type <- old_set_prec_type session
      | LSI.Pcg ->
          c_pcg session maxl prec_type;
          compat.set_maxl <- old_set_maxl session;
          compat.set_prec_type <- old_set_prec_type session
      | _ -> assert false

    let solver (type s)
          (LSI.(LS ({ rawptr; solver; compat } as hls) as ls))
          ?jac_times_vec ?jac_times_rhs (prec_type, set_prec) session nv =
      let jac_times_setup, jac_times_vec =
        match jac_times_vec with
        | None -> None, None
        | Some _ when jac_times_rhs <> None ->
            invalid_arg "cannot pass both jac_times_vec and jac_times_rhs"
        | Some (ojts, jtv) -> ojts, Some jtv
      in
      if Sundials_impl.Version.lt530 && jac_times_rhs <> None
        then raise Config.NotImplementedBySundialsVersion;
      if Sundials_impl.Version.in_compat_mode2 then begin
        if jac_times_setup <> None then
          raise Config.NotImplementedBySundialsVersion;
        make_compat compat prec_type solver session;
        session.ls_solver <- LSI.HLS hls;
        set_prec session nv;
        session.ls_callbacks <- SpilsCallback1 (jac_times_vec, None);
        if jac_times_vec <> None then c_set_jac_times session true false
      end else
        if Sundials_impl.Version.in_compat_mode2_3 then c_spils_set_linear_solver session rawptr
        else c_set_linear_solver session rawptr None false false;
        LSI.attach ls;
        session.ls_solver <- LSI.HLS hls;
        LSI.(impl_set_prec_type rawptr solver prec_type false);
        set_prec session nv;
        match jac_times_rhs with
        | Some jtrhsfn -> begin
            session.ls_callbacks <- SpilsCallback2 jtrhsfn;
            c_set_jac_times_rhsfn session true
          end
        | None -> begin
            session.ls_callbacks <-
              SpilsCallback1 (jac_times_vec, jac_times_setup);
            if jac_times_setup <> None || jac_times_vec <> None then
              c_set_jac_times session (jac_times_setup <> None)
                                      (jac_times_vec <> None)
          end

    let set_jac_times s ?jac_times_setup f =
      if Sundials_impl.Version.in_compat_mode2 && jac_times_setup <> None then
          raise Config.NotImplementedBySundialsVersion;
      match s.ls_callbacks with
      | SpilsCallback1 _ ->
          c_set_jac_times s (jac_times_setup <> None) true;
          s.ls_callbacks <- SpilsCallback1 (Some f, jac_times_setup)
      | _ -> raise LinearSolver.InvalidLinearSolver

    let clear_jac_times s =
      match s.ls_callbacks with
      | SpilsCallback1 _ ->
          c_set_jac_times s false false;
          s.ls_callbacks <- SpilsCallback1 (None, None)
      | _ -> raise LinearSolver.InvalidLinearSolver

    let set_preconditioner s ?setup solve =
      match s.ls_callbacks with
      | SpilsCallback1 _ | SpilsCallback2 _ ->
          c_set_preconditioner s (setup <> None);
          s.ls_precfns <- PrecFns { prec_setup_fn = setup;
                                    prec_solve_fn = solve }
      | _ -> raise LinearSolver.InvalidLinearSolver

    external set_jac_eval_frequency : ('a, 'k) session -> int -> unit
        = "sunml_arkode_ark_set_jac_eval_frequency"

    let set_jac_eval_frequency s maxsteps =
      if Sundials_impl.Version.in_compat_mode2_3 then ls_check_spils s;
      set_jac_eval_frequency s maxsteps

    external set_linear_solution_scaling : ('d, 'k) session -> bool -> unit
      = "sunml_arkode_ark_set_linear_solution_scaling"

    external set_eps_lin            : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_eps_lin"

    let set_eps_lin s epsl =
      if Sundials_impl.Version.in_compat_mode2_3 then ls_check_spils s;
      set_eps_lin s epsl

    external set_ls_norm_factor : ('d, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_ls_norm_factor"

    external get_num_lin_iters      : ('a, 'k) session -> int
      = "sunml_arkode_ark_get_num_lin_iters"

    let get_num_lin_iters s =
      if Sundials_impl.Version.in_compat_mode2_3 then ls_check_spils s;
      get_num_lin_iters s

    external get_num_lin_conv_fails : ('a, 'k) session -> int
      = "sunml_arkode_ark_get_num_lin_conv_fails"

    let get_num_lin_conv_fails s =
      if Sundials_impl.Version.in_compat_mode2_3 then ls_check_spils s;
      get_num_lin_conv_fails s

    external get_work_space         : ('a, 'k) session -> int * int
      = "sunml_arkode_spils_get_work_space"

    let get_work_space s =
      if Sundials_impl.Version.in_compat_mode2_3 then ls_check_spils s;
      get_work_space s

    external get_num_prec_evals     : ('a, 'k) session -> int
      = "sunml_arkode_ark_get_num_prec_evals"

    let get_num_prec_evals s =
      if Sundials_impl.Version.in_compat_mode2_3 then ls_check_spils s;
      get_num_prec_evals s

    external get_num_prec_solves    : ('a, 'k) session -> int
      = "sunml_arkode_ark_get_num_prec_solves"

    let get_num_prec_solves s =
      ls_check_spils s;
      get_num_prec_solves s

    external get_num_jtsetup_evals   : ('a, 'k) session -> int
      = "sunml_arkode_ark_get_num_jtsetup_evals"

    let get_num_jtsetup_evals s =
      ls_check_spils s;
      get_num_jtsetup_evals s

    external get_num_jtimes_evals   : ('a, 'k) session -> int
      = "sunml_arkode_ark_get_num_jtimes_evals"

    let get_num_jtimes_evals s =
      ls_check_spils s;
      get_num_jtimes_evals s

    external get_num_lin_rhs_evals  : ('a, 'k) session -> int
      = "sunml_arkode_spils_get_num_lin_rhs_evals"

    let get_num_lin_rhs_evals s =
      ls_check_spils s;
      get_num_lin_rhs_evals s

  end (* }}} *)

  module Mass = struct (* {{{ *)
    include MassTypes

    module Dls = struct (* {{{ *)
      include MassTypes.Direct'
      include LinearSolver.Direct

      (* Sundials < 3.0.0 *)
      external c_dls_mass_dense : 'k serial_session -> int -> unit
        = "sunml_arkode_dls_mass_dense"

      (* Sundials < 3.0.0 *)
      external c_dls_mass_lapack_dense
        : 'k serial_session -> int -> unit
        = "sunml_arkode_dls_mass_lapack_dense"

      (* Sundials < 3.0.0 *)
      external c_dls_mass_band
        : 'k serial_session -> int -> int -> int -> unit
        = "sunml_arkode_dls_mass_band"

      (* Sundials < 3.0.0 *)
      external c_dls_mass_lapack_band
        : 'k serial_session -> int -> int -> int -> unit
        = "sunml_arkode_dls_mass_lapack_band"

      (* Sundials < 3.0.0 *)
      external c_mass_klu
        : 'k serial_session -> 's Matrix.Sparse.sformat
          -> int -> int -> unit
        = "sunml_arkode_mass_klu_init"

      (* Sundials < 3.0.0 *)
      external c_klu_set_ordering
        : 'k serial_session -> LSI.Klu.ordering -> unit
        = "sunml_arkode_mass_klu_set_ordering"

      (* Sundials < 3.0.0 *)
      external c_klu_reinit
        : 'k serial_session -> int -> int -> unit
        = "sunml_arkode_mass_klu_reinit"

      (* Sundials < 3.0.0 *)
      external c_mass_superlumt
        : 'k serial_session -> int -> int -> int -> unit
        = "sunml_arkode_mass_superlumt_init"

      (* Sundials < 3.0.0 *)
      external c_superlumt_set_ordering
        : 'k serial_session -> LSI.Superlumt.ordering -> unit
        = "sunml_arkode_mass_superlumt_set_ordering"

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
            (solver_data : (s, 'nd, 'nk, tag) LSI.solver_data)
            (mat : ('k, s, 'nd, 'nk) Matrix.t) session =
        match solver_data with
        | LSI.Dense ->
            let m, n = Matrix.(Dense.size (unwrap mat)) in
            if m <> n then raise LinearSolver.MatrixNotSquare;
            c_dls_mass_dense session m
        | LSI.LapackDense ->
            let m, n = Matrix.(Dense.size (unwrap mat)) in
            if m <> n then raise LinearSolver.MatrixNotSquare;
            c_dls_mass_lapack_dense session m

        | LSI.Band ->
            let open Matrix.Band in
            let { n; mu; ml } = dims (Matrix.unwrap mat) in
            c_dls_mass_band session n mu ml
        | LSI.LapackBand ->
            let open Matrix.Band in
            let { n; mu; ml } = dims (Matrix.unwrap mat) in
            c_dls_mass_lapack_band session n mu ml

        | LSI.Klu sinfo ->
            if not Config.klu_enabled
              then raise Config.NotImplementedBySundialsVersion;
            let smat = Matrix.unwrap mat in
            let m, n = Matrix.Sparse.size smat in
            let nnz, _ = Matrix.Sparse.dims smat in
            if m <> n then raise LinearSolver.MatrixNotSquare;
            let open LSI.Klu in
            sinfo.set_ordering <- klu_set_ordering session;
            sinfo.reinit <- klu_reinit session;
            c_mass_klu session (Matrix.Sparse.sformat smat) m nnz;
            (match sinfo.ordering with None -> ()
                                     | Some o -> c_klu_set_ordering session o)

        | LSI.Superlumt sinfo ->
            if not Config.superlumt_enabled
              then raise Config.NotImplementedBySundialsVersion;
            let smat = Matrix.unwrap mat in
            let m, n = Matrix.Sparse.size smat in
            let nnz, _ = Matrix.Sparse.dims smat in
            if m <> n then raise LinearSolver.MatrixNotSquare;
            let open LSI.Superlumt in
            sinfo.set_ordering <- superlumt_set_ordering session;
            c_mass_superlumt session m nnz sinfo.num_threads;
            (match sinfo.ordering with None -> ()
                                     | Some o -> c_superlumt_set_ordering session o)

      | _ -> assert false

      let set_mass_callbacks (type m) (type tag)
            (massfn : m mass_fn)
            (solver_data : (m, 'nd, 'nk, tag) LSI.solver_data)
            (mat : ('mk, m, 'nd, 'nk) Matrix.t) session =
        let cb = { massfn = massfn; mmat  = (None : m option) } in
        begin match solver_data with
        | LSI.Dense ->
            session.mass_callbacks <-
              DlsDenseMassCallback (cb, Matrix.unwrap mat)
        | LSI.LapackDense ->
            session.mass_callbacks <-
              DlsDenseMassCallback (cb, Matrix.unwrap mat)
        | LSI.Band ->
            session.mass_callbacks <-
              DlsBandMassCallback (cb, Matrix.unwrap mat)
        | LSI.LapackBand ->
            session.mass_callbacks <-
              DlsBandMassCallback (cb, Matrix.unwrap mat)
        | LSI.Klu _ ->
            session.mass_callbacks <-
              SlsKluMassCallback (cb, Matrix.unwrap mat)
        | LSI.Superlumt _ ->
            session.mass_callbacks
              <- SlsSuperlumtMassCallback (cb, Matrix.unwrap mat)
        | LSI.Custom _ ->
            session.mass_callbacks
              <- DirectCustomMassCallback (cb, Matrix.unwrap mat)
        | _ -> assert false
        end;
        session.mass_precfns <- NoMassPrecFns

      (* 3.0.0 <= Sundials < 4.0.0 *)
      external c_dls_set_mass_linear_solver
        : 'k serial_session
        -> ('m, Nvector_serial.data, 'k) LSI.cptr
          -> ('mk, 'm, Nvector_serial.data, 'k) Matrix.t
          -> bool
          -> unit
        = "sunml_arkode_dls_set_mass_linear_solver"

      (* 4.0.0 <= Sundials *)
      external c_set_mass_linear_solver
        : ('d, 'k) session
          -> ('m, 'd, 'k) LSI.cptr
          -> ('mk, 'm, 'd, 'k) Matrix.t option
          -> bool
          -> unit
        = "sunml_arkode_ark_set_mass_linear_solver"

      (* 4.0.0 <= Sundials *)
      external c_set_mass_fn : ('a, 'k) session -> unit
        = "sunml_arkode_ark_set_mass_fn"

      let assert_matrix = function
        | Some m -> m
        | None -> failwith "a direct linear solver is required"

      let solver massfn time_dep
                 (LSI.(LS ({ rawptr; solver; matrix } as hls)) as ls)
                 session nv =
        let matrix = assert_matrix matrix in
        set_mass_callbacks massfn solver matrix session;
        if Sundials_impl.Version.in_compat_mode2 then make_compat solver matrix session
        else if Sundials_impl.Version.in_compat_mode2_3
             then c_dls_set_mass_linear_solver session rawptr matrix time_dep
        else (c_set_mass_linear_solver session rawptr (Some matrix) time_dep;
              c_set_mass_fn session);
        LSI.attach ls;
        session.mass_solver <- LSI.HLS hls

      (* Sundials < 3.0.0 *)
      let invalidate_callback session =
        if Sundials_impl.Version.in_compat_mode2 then
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
          = "sunml_arkode_dls_get_mass_work_space"

      let get_work_space s =
        mass_check_direct s;
        get_work_space s

      external c_get_num_mass_setups : 'k serial_session -> int
        = "sunml_arkode_ark_get_num_mass_setups"

      let get_num_setups s =
        mass_check_direct s;
        c_get_num_mass_setups s

      external get_num_mult_setups : 'k serial_session -> int
        = "sunml_arkode_ark_get_num_mass_mult_setups"

      external c_get_num_mass_solves
        : 'k serial_session -> int
        = "sunml_arkode_ark_get_num_mass_solves"

      let get_num_solves s =
        mass_check_direct s;
        c_get_num_mass_solves s

      external c_get_num_mass_mult
        : 'k serial_session -> int
        = "sunml_arkode_dls_get_num_mass_mult"

      let get_num_mult s =
        mass_check_direct s;
        c_get_num_mass_mult s

    end (* }}} *)

    module Spils = struct (* {{{ *)
      include MassTypes.Iterative'
      include LinearSolver.Iterative

      (* Sundials < 3.0.0 *)
      external c_spgmr
        : ('a, 'k) session
          -> int -> LSI.Iterative.preconditioning_type -> unit
        = "sunml_arkode_spils_mass_spgmr"

      (* Sundials < 3.0.0 *)
      external c_spbcgs
        : ('a, 'k) session
          -> int -> LSI.Iterative.preconditioning_type -> unit
        = "sunml_arkode_spils_mass_spbcgs"

      (* Sundials < 3.0.0 *)
      external c_sptfqmr
        : ('a, 'k) session
          -> int -> LSI.Iterative.preconditioning_type -> unit
        = "sunml_arkode_spils_mass_sptfqmr"

      (* Sundials < 3.0.0 *)
      external c_spfgmr
        : ('a, 'k) session
          -> int -> LSI.Iterative.preconditioning_type -> unit
        = "sunml_arkode_spils_mass_spfgmr"

      (* Sundials < 3.0.0 *)
      external c_pcg
        : ('a, 'k) session
          -> int -> LSI.Iterative.preconditioning_type -> unit
        = "sunml_arkode_spils_mass_pcg"

      (* Sundials < 3.0.0 *)
      external c_set_gs_type
        : ('a, 'k) session
          -> LSI.Iterative.gramschmidt_type -> unit
        = "sunml_arkode_spils_set_mass_gs_type"

      (* Sundials < 3.0.0 *)
      external c_set_maxl : ('a, 'k) session -> int -> unit
        = "sunml_arkode_spils_set_mass_maxl"

      (* Sundials < 3.0.0 *)
      external c_set_prec_type
        : ('a, 'k) session -> LSI.Iterative.preconditioning_type -> unit
        = "sunml_arkode_spils_set_mass_prec_type"

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
        = "sunml_arkode_ark_set_mass_times"

      external c_set_preconditioner
        : ('a, 'k) session -> bool -> unit
        = "sunml_arkode_ark_set_mass_preconditioner"

      (* Sundials < 4.0.0 *)
      external c_spils_set_mass_linear_solver
        : ('a, 'k) session
          -> ('m, 'a, 'k) LSI.cptr
          -> bool
          -> bool
          -> unit
        = "sunml_arkode_spils_set_mass_linear_solver"

      (* 4.0.0 <= Sundials *)
      external c_set_mass_linear_solver
        : ('d, 'k) session
          -> ('m, 'd, 'k) LSI.cptr
          -> ('mk, 'm, 'd, 'k) Matrix.t option
          -> bool
          -> unit
        = "sunml_arkode_ark_set_mass_linear_solver"

      let init_preconditioner solve setup session nv =
        c_set_preconditioner session (setup <> None);
        session.mass_precfns <- MassPrecFns { prec_solve_fn = solve;
                                              prec_setup_fn = setup }

      let prec_none = LSI.Iterative.(PrecNone,
                        fun session nv -> session.mass_precfns <- NoMassPrecFns)

      let prec_left ?setup solve  = LSI.Iterative.(PrecLeft,
                                                init_preconditioner solve setup)

      let prec_right ?setup solve = LSI.Iterative.(PrecRight,
                                                init_preconditioner solve setup)

      let prec_both ?setup solve  = LSI.Iterative.(PrecBoth,
                                                init_preconditioner solve setup)

      (* Sundials < 3.0.0 *)
      let make_compat (type tag)
            (LSI.Iterative.({ maxl; gs_type }) as compat)
            prec_type
            (solver_data : ('s, 'nd, 'nk, tag) LSI.solver_data) session =
        (match solver_data with
         | LSI.Spgmr ->
             c_spgmr session maxl prec_type;
             (match gs_type with None -> () | Some t -> c_set_gs_type session t);
             compat.set_gs_type <- old_set_gs_type session;
             compat.set_prec_type <- old_set_prec_type session
         | LSI.Spfgmr ->
             c_spfgmr session maxl prec_type;
             (match gs_type with None -> () | Some t -> c_set_gs_type session t);
             compat.set_gs_type <- old_set_gs_type session;
             compat.set_prec_type <- old_set_prec_type session
         | LSI.Spbcgs ->
             c_spbcgs session maxl prec_type;
             compat.set_maxl <- old_set_maxl session;
             compat.set_prec_type <- old_set_prec_type session
         | LSI.Sptfqmr ->
             c_sptfqmr session maxl prec_type;
             compat.set_maxl <- old_set_maxl session;
             compat.set_prec_type <- old_set_prec_type session
         | LSI.Pcg ->
             c_pcg session maxl prec_type;
             compat.set_maxl <- old_set_maxl session;
             compat.set_prec_type <- old_set_prec_type session
         | _ -> assert false)

      let solver (type s)
            LSI.(LS { rawptr; solver; compat })
            ?mass_times_setup mass_times_vec time_dep (prec_type, set_prec)
            session nv =
        if Sundials_impl.Version.in_compat_mode2 then begin
          if mass_times_setup <> None then
            raise Config.NotImplementedBySundialsVersion;
          make_compat compat prec_type solver session;
          set_prec session nv;
          session.mass_callbacks <- SpilsMassCallback (mass_times_vec, None)
        end else
          if Sundials_impl.Version.in_compat_mode2_3
             then c_spils_set_mass_linear_solver session rawptr
                                                 (mass_times_setup <> None)
                                                 time_dep
          else (c_set_mass_linear_solver session rawptr None false;
                c_set_mass_times session (mass_times_setup <> None));
          LSI.(impl_set_prec_type rawptr solver prec_type false);
          set_prec session nv;
          session.mass_callbacks <- SpilsMassCallback
                                      (mass_times_vec, mass_times_setup)

      let set_times s ?mass_times_setup mass_times_vec =
        if Sundials_impl.Version.in_compat_mode2 && mass_times_setup <> None then
            raise Config.NotImplementedBySundialsVersion;
        match s.mass_callbacks with
        | SpilsMassCallback _ ->
            c_set_mass_times s (mass_times_setup <> None);
            s.mass_callbacks <- SpilsMassCallback
                                  (mass_times_vec, mass_times_setup)
        | _ -> raise LinearSolver.InvalidLinearSolver

      let set_preconditioner s ?setup solve =
        match s.mass_callbacks with
        | SpilsMassCallback _ ->
            c_set_preconditioner s (setup <> None);
            s.mass_precfns <- MassPrecFns { prec_setup_fn = setup;
                                            prec_solve_fn = solve }
        | _ -> raise LinearSolver.InvalidLinearSolver

      external set_eps_lin            : ('a, 'k) session -> float -> unit
          = "sunml_arkode_ark_set_mass_eps_lin"

      let set_eps_lin s epsl =
        mass_check_spils s;
        set_eps_lin s epsl

      external set_ls_norm_factor : ('d, 'k) session -> float -> unit
        = "sunml_arkode_ark_set_mass_ls_norm_factor"

      external get_num_lin_iters      : ('a, 'k) session -> int
          = "sunml_arkode_ark_get_num_mass_iters"

      let get_num_lin_iters s =
        mass_check_spils s;
        get_num_lin_iters s

      external get_num_conv_fails     : ('a, 'k) session -> int
          = "sunml_arkode_ark_get_num_mass_conv_fails"

      let get_num_conv_fails s =
        mass_check_spils s;
        get_num_conv_fails s

      external get_num_mtsetups   : ('a, 'k) session -> int
        = "sunml_arkode_ark_get_num_mtsetups"

      let get_num_mtsetups s =
        mass_check_spils s;
        get_num_mtsetups s

      external get_num_mass_mult       : ('a, 'k) session -> int
          = "sunml_arkode_spils_get_num_mass_mult"

      let get_num_mass_mult s =
        mass_check_spils s;
        get_num_mass_mult s

      external get_work_space         : ('a, 'k) session -> int * int
          = "sunml_arkode_spils_get_mass_work_space"

      let get_work_space s =
        mass_check_spils s;
        get_work_space s

      external get_num_prec_evals     : ('a, 'k) session -> int
          = "sunml_arkode_ark_get_num_mass_prec_evals"

      let get_num_prec_evals s =
        mass_check_spils s;
        get_num_prec_evals s

      external get_num_prec_solves    : ('a, 'k) session -> int
          = "sunml_arkode_ark_get_num_mass_prec_solves"

      let get_num_prec_solves s =
        mass_check_spils s;
        get_num_prec_solves s
    end (* }}} *)
  end (* }}} *)

  external sv_tolerances
      : ('a, 'k) session -> float -> ('a, 'k) nvector -> unit
      = "sunml_arkode_ark_sv_tolerances"
  external ss_tolerances
      : ('a, 'k) session -> float -> float -> unit
      = "sunml_arkode_ark_ss_tolerances"
  external wf_tolerances 
      : ('a, 'k) session -> unit
      = "sunml_arkode_ark_wf_tolerances"

  let set_tolerances s tol =
    match tol with
    | SStolerances (rel, abs) -> (s.errw <- dummy_errw; ss_tolerances s rel abs)
    | SVtolerances (rel, abs) -> (if Sundials_configuration.safe then s.checkvec abs;
                                  s.errw <- dummy_errw; sv_tolerances s rel abs)
    | WFtolerances ferrw -> (s.errw <- ferrw; wf_tolerances s)

  external ress_tolerance  : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_ress_tolerance"
  external resv_tolerance  : ('a, 'k) session -> ('a, 'k) nvector -> unit
      = "sunml_arkode_ark_resv_tolerance"
  external resf_tolerance  : ('a, 'k) session -> unit
      = "sunml_arkode_ark_resf_tolerance"

  type 'data res_weight_fun = 'data -> 'data -> unit

  type ('data, 'kind) res_tolerance =
    | ResStolerance of float
    | ResVtolerance of ('data, 'kind) Nvector.t
    | ResFtolerance of 'data res_weight_fun

  let set_res_tolerance s tol =
    match tol with
    | ResStolerance abs -> (s.resw <- dummy_resw;
                            s.uses_resv <- false;
                            ress_tolerance s abs)
    | ResVtolerance abs -> (if Sundials_configuration.safe then s.checkvec abs;
                            s.uses_resv <- true;
                            s.resw <- dummy_resw;
                            resv_tolerance s abs)
    | ResFtolerance fresw -> (s.resw <- fresw;
                              s.uses_resv <- false;
                              resf_tolerance s)

  (* Sundials < 4.0.0 *)
  external set_fixed_point : ('a, 'k) session -> int -> unit
    = "sunml_arkode_ark_set_fixed_point"

  (* Sundials < 4.0.0 *)
  external set_newton : ('a, 'k) session -> unit
    = "sunml_arkode_ark_set_newton"

  external set_linear : ('a, 'k) session -> bool -> unit
    = "sunml_arkode_ark_set_linear"

  external set_nonlinear : ('a, 'k) session -> unit
    = "sunml_arkode_ark_set_nonlinear"

  external c_set_order : ('a, 'k) session -> int -> unit
    = "sunml_arkode_ark_set_order"

  external c_session_finalize : ('a, 'k) session -> unit
      = "sunml_arkode_ark_session_finalize"

  let session_finalize s =
    Dls.invalidate_callback s;
    (match s.nls_solver with
     | None -> ()
     | Some nls -> NLSI.detach nls);
    c_session_finalize s

  (* 4.0.0 <= Sundials *)
  external c_set_nonlinear_solver
      : ('d, 'k) session
        -> ('d, 'k, ('d, 'k) session, [`Nvec]) NLSI.cptr
        -> unit
      = "sunml_arkode_ark_set_nonlinear_solver"

  external c_init :
    ('a, 'k) session Weak.t
    -> bool             (* f_i given *)
    -> bool             (* f_e given *)
    -> ('a, 'k) nvector (* y_0 *)
    -> float            (* t_0 *)
    -> (arkstep arkode_mem * c_weak_ref)
    = "sunml_arkode_ark_init"

  let init prob tol ?restol ?order ?mass ?(roots=no_roots) t0 y0 =
    let (nroots, roots) = roots in
    let checkvec = Nvector.check y0 in
    if Sundials_configuration.safe && nroots < 0
      then invalid_arg "number of root functions is negative";
    let problem, fi, fe, nlsolver, lsolver, lin =
      match prob with
      | Implicit { irhsfn=fi; linearity=l; nlsolver=nls; lsolver=ls } ->
          ImplicitOnly,        Some fi, None,    nls,  ls,   Some l
      | Explicit fe ->
          ExplicitOnly,        None,    Some fe, None, None, None
      | ImEx (fe, { irhsfn=fi; linearity=l; nlsolver=nls; lsolver=ls }) ->
          ImplicitAndExplicit, Some fi, Some fe, nls,  ls,   Some l
    in
    let weakref = Weak.create 1 in
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
            rhsfn1       = (match fi with Some f -> f | None -> dummy_rhsfn1);
            rhsfn2       = (match fe with Some f -> f | None -> dummy_rhsfn2);

            rootsfn      = roots;
            errh         = dummy_errh;
            errw         = dummy_errw;
            resw         = dummy_resw;

            adaptfn      = dummy_adaptfn;
            stabfn       = dummy_stabfn;
            resizefn     = dummy_resizefn;
            poststepfn   = dummy_poststepfn;
            stagepredictfn = dummy_stagepredictfn;
            preinnerfn   = dummy_preinnerfn;
            postinnerfn  = dummy_postinnerfn;
            preinnerarray = empty_preinnerarray ();

            linsolver      = None;
            ls_solver      = LSI.NoHLS;
            ls_callbacks   = NoCallbacks;
            ls_precfns     = NoPrecFns;

            mass_solver    = LSI.NoHLS;
            mass_callbacks = NoMassCallbacks;
            mass_precfns   = NoMassPrecFns;

            nls_solver     = None;

            inner_session  = None;
          } in
    Gc.finalise session_finalize session;
    Weak.set weakref 0 (Some session);
    (* Now the session is safe to use.  If any of the following fails and raises
       an exception, the GC will take care of freeing arkode_mem and backref.  *)
    if nroots > 0 then
      c_root_init session nroots;
    set_tolerances session tol;
    if Sundials_impl.Version.in_compat_mode2_3 then begin
      match nlsolver with
      | Some NLSI.{ solver = NewtonSolver _ } | None -> () (* the default *)
      | Some NLSI.{ solver = FixedPointSolver (_, fpm) } ->
          set_fixed_point session fpm
      | _ -> raise Config.NotImplementedBySundialsVersion
    end else begin
      match nlsolver with
      | Some ({ NLSI.rawptr = nlcptr } as nls) ->
          NLSI.attach nls;
          session.nls_solver <- Some nls;
          c_set_nonlinear_solver session nlcptr
      | _ -> ()
    end;
    (match lsolver with
     | None -> ()
     | Some ls -> session.linsolver <- Some ls;
                  ls session y0);
    (match lin with
     | Some (Linear timedepend) -> set_linear session timedepend
     | _ -> ());
    (match restol with Some rtol -> set_res_tolerance session rtol | None -> ());
    (match order with Some o -> c_set_order session o | None -> ());
    (match mass with Some msolver -> msolver session y0 | None -> ());
    session

  let get_num_roots { nroots } = nroots

  external reset : ('d, 'k) session -> float -> ('d, 'k) Nvector.t
      = "sunml_arkode_ark_reset"

  external c_reinit
      : ('a, 'k) session -> float -> ('a, 'k) nvector -> unit
      = "sunml_arkode_ark_reinit"

  let reinit session ?problem ?order ?mass ?roots t0 y0 =
    if Sundials_configuration.safe then session.checkvec y0;
    Dls.invalidate_callback session;
    let lin, nlsolver, lsolver =
      match problem with
      | None -> None, None, None
      | Some (Implicit { irhsfn = fi; linearity; nlsolver; lsolver }) ->
          session.problem <- ImplicitOnly;
          session.rhsfn1 <- fi;
          session.rhsfn2 <- dummy_rhsfn2;
          Some linearity, nlsolver, lsolver
      | Some (Explicit fe) ->
          session.problem <- ExplicitOnly;
          session.rhsfn1 <- dummy_rhsfn1;
          session.rhsfn2 <- fe;
          None, None, None
      | Some (ImEx (fe, { irhsfn = fi; linearity; nlsolver; lsolver })) ->
          session.problem <- ImplicitAndExplicit;
          session.rhsfn1 <- fi;
          session.rhsfn2 <- fe;
          Some linearity, nlsolver, lsolver
    in
    (match lin with
     | None -> ()
     | Some Nonlinear -> set_nonlinear session
     | Some (Linear timedepend) -> set_linear session timedepend);
    (if Sundials_impl.Version.in_compat_mode2_3 then begin
      match nlsolver with
      | None -> ()
      | Some NLSI.{ solver = NewtonSolver _ } -> set_newton session
      | Some NLSI.{ solver = FixedPointSolver (_, fpm) } ->
          set_fixed_point session fpm
      | _ -> raise Config.NotImplementedBySundialsVersion
    end else begin
      match nlsolver with
      | Some ({ NLSI.rawptr = nlcptr } as nls) ->
          (match session.nls_solver with
           | None -> () | Some old_nls -> NLSI.detach old_nls);
          NLSI.attach nls;
          session.nls_solver <- Some nls;
          c_set_nonlinear_solver session nlcptr
      | _ -> ()
    end);
    (match lsolver with None -> () | Some ls -> ls session y0);
    session.linsolver <- lsolver;
    (match mass with Some msolver -> msolver session y0 | None -> ());
    c_reinit session t0 y0;
    (match order with Some o -> c_set_order session o | None -> ());
    (match roots with
     | None -> ()
     | Some roots -> root_init session roots)

  external c_resize
      : ('a, 'k) session -> bool -> float -> float -> ('a, 'k) nvector -> unit
      = "sunml_arkode_ark_resize"

  let resize session ?resize_nvec ?lsolver ?mass tol ?restol hscale ynew t0 =
    session.checkvec <- Nvector.check ynew;
    (match lsolver with None -> () | ls -> session.linsolver <- ls);
    (match resize_nvec with None -> () | Some f -> session.resizefn <- f);
    c_resize session (resize_nvec <> None) hscale t0 ynew;
    session.resizefn <- dummy_resizefn;
    (match Config.sundials_version with
     | 2,6,1 | 2,6,2 -> () (* avoid a segmentation fault in earlier versions *)
     | _ -> set_tolerances session tol);
    (match restol with
     | Some rt -> set_res_tolerance session rt
     | _ when session.uses_resv -> set_res_tolerance session (ResStolerance 1.e-9)
     | _ -> ());
    (match session.linsolver with Some ls -> ls session ynew | None -> ());
    (match mass with Some msolver -> msolver session ynew | None -> ())

  external get_root_info  : ('a, 'k) session -> Roots.t -> unit
      = "sunml_arkode_ark_get_root_info"

  external c_solve_normal : ('a, 'k) session -> float -> ('a, 'k) nvector
                                -> float * solver_result
      = "sunml_arkode_ark_solve_normal"

  let solve_normal s t y =
    if Sundials_configuration.safe then s.checkvec y;
    c_solve_normal s t y

  external c_solve_one_step : ('a, 'k) session -> float -> ('a, 'k) nvector
                                -> float * solver_result
      = "sunml_arkode_ark_solve_one_step"

  let solve_one_step s t y =
    if Sundials_configuration.safe then s.checkvec y;
    c_solve_one_step s t y

  external c_get_dky
      : ('a, 'k) session -> float -> int -> ('a, 'k) nvector -> unit
      = "sunml_arkode_ark_get_dky"

  let get_dky s y =
    if Sundials_configuration.safe then s.checkvec y;
    fun t k -> c_get_dky s t k y

  (* Synchronized with arkode_timestepper_stats_index in arkode_ml.h *)
  type timestepper_stats = {
      exp_steps           : int;
      acc_steps           : int;
      step_attempts       : int;
      num_nfe_evals       : int;
      num_nfi_evals       : int;
      num_lin_solv_setups : int;
      num_err_test_fails  : int;
    }

  external get_timestepper_stats : ('a, 'k) session -> timestepper_stats
      = "sunml_arkode_ark_get_timestepper_stats"

  external get_step_stats : ('a, 'k) session -> step_stats
      = "sunml_arkode_ark_get_step_stats"

  external get_work_space         : ('a, 'k) session -> int * int
      = "sunml_arkode_ark_get_work_space"

  external get_num_steps          : ('a, 'k) session -> int
      = "sunml_arkode_ark_get_num_steps"

  external get_num_exp_steps      : ('d, 'k) session -> int
      = "sunml_arkode_ark_get_num_exp_steps"

  external get_num_acc_steps      : ('d, 'k) session -> int
      = "sunml_arkode_ark_get_num_acc_steps"

  external get_num_step_attempts  : ('d, 'k) session -> int
      = "sunml_arkode_ark_get_num_step_attempts"

  external get_num_rhs_evals      : ('a, 'k) session -> int * int
      = "sunml_arkode_ark_get_num_rhs_evals"

  external get_num_lin_solv_setups : ('a, 'k) session -> int
      = "sunml_arkode_ark_get_num_lin_solv_setups"

  external get_num_err_test_fails : ('a, 'k) session -> int
      = "sunml_arkode_ark_get_num_err_test_fails"

  external get_actual_init_step   : ('a, 'k) session -> float
      = "sunml_arkode_ark_get_actual_init_step"

  external get_last_step          : ('a, 'k) session -> float
      = "sunml_arkode_ark_get_last_step"

  external get_current_step       : ('a, 'k) session -> float
      = "sunml_arkode_ark_get_current_step"

  external get_current_time       : ('a, 'k) session -> float
      = "sunml_arkode_ark_get_current_time"

  external get_current_state : ('d, 'k) session -> 'd
      = "sunml_arkode_ark_get_current_state"

  external get_nonlin_system_data
    : ('d, 'k) session -> 'd nonlin_system_data
    = "sunml_arkode_ark_get_nonlin_system_data"

  external compute_state
    : ('d, 'k) session
      -> ('d, 'k) Nvector.t
      -> ('d, 'k) Nvector.t
      -> unit
    = "sunml_arkode_ark_compute_state"

  let compute state s ycor yn =
    if Sundials_configuration.safe then (s.checkvec ycor; s.checkvec yn);
    compute_state s ycor yn

  external get_current_gamma : ('d, 'k) session -> float
      = "sunml_arkode_ark_get_current_gamma"

  let print_timestepper_stats s oc =
    let stats = get_timestepper_stats s
    in
    Printf.fprintf oc "exp_steps = %d\n"           stats.exp_steps;
    Printf.fprintf oc "acc_steps = %d\n"           stats.acc_steps;
    Printf.fprintf oc "step_attempts = %d\n"       stats.step_attempts;
    Printf.fprintf oc "num_nfe_evals = %d\n"       stats.num_nfe_evals;
    Printf.fprintf oc "num_nfi_evals = %d\n"       stats.num_nfi_evals;
    Printf.fprintf oc "num_lin_solv_setups = %d\n" stats.num_lin_solv_setups;
    Printf.fprintf oc "num_err_test_fails = %d\n"  stats.num_err_test_fails

  let print_step_stats s oc =
    print_step_stats oc (get_step_stats s)

  external set_diagnostics : ('a, 'k) session -> Logfile.t -> unit
      = "sunml_arkode_ark_set_diagnostics"

  external clear_diagnostics : ('a, 'k) session -> unit
      = "sunml_arkode_ark_clear_diagnostics"

  external set_error_file : ('a, 'k) session -> Logfile.t -> unit
      = "sunml_arkode_ark_set_error_file"

  external c_set_err_handler_fn  : ('a, 'k) session -> unit
      = "sunml_arkode_ark_set_err_handler_fn"

  let set_err_handler_fn s ferrh =
    s.errh <- ferrh;
    c_set_err_handler_fn s

  external clear_err_handler_fn  : ('a, 'k) session -> unit
      = "sunml_arkode_ark_clear_err_handler_fn"

  let clear_err_handler_fn s =
    s.errh <- dummy_errh;
    clear_err_handler_fn s

  external c_set_imex             : ('a, 'k) session -> unit
      = "sunml_arkode_ark_set_imex"

  external c_set_explicit         : ('a, 'k) session -> unit
      = "sunml_arkode_ark_set_explicit"

  external c_set_implicit         : ('a, 'k) session -> unit
      = "sunml_arkode_ark_set_implicit"

  let set_imex s =
    (if s.rhsfn1 == dummy_rhsfn1 || s.rhsfn2 == dummy_rhsfn2 then raise IllInput);
    s.problem <- ImplicitAndExplicit;
    c_set_imex s

  let set_explicit s =
    (if s.rhsfn2 == dummy_rhsfn2 then raise IllInput);
    s.problem <- ExplicitOnly;
    c_set_explicit s

  let set_implicit s =
    (if s.rhsfn1 == dummy_rhsfn1 then raise IllInput);
    s.problem <- ImplicitOnly;
    c_set_implicit s

  external c_set_tables
    : ('d, 'k) session
      -> int (* q *)
      -> int (* p *)
      -> ButcherTable.t option
      -> ButcherTable.t option
      -> unit
      = "sunml_arkode_ark_set_tables"

  let set_tables s ?(global_method_order=(-1))
                   ?(global_embedding_order=(-1))
                   ?implicit_table ?explicit_table () =
    if (implicit_table <> None && s.rhsfn1 == dummy_rhsfn1)
       || (explicit_table <> None && s.rhsfn2 == dummy_rhsfn2)
    then raise IllInput;
    if implicit_table = None && explicit_table = None then raise IllInput;
    if not Sundials_impl.Version.in_compat_mode2_3
       && implicit_table <> None
       && explicit_table <> None
       && (global_method_order < 0 || global_embedding_order < 0)
    then raise IllInput;
    c_set_tables s global_method_order global_embedding_order
                   implicit_table explicit_table

  external c_set_table_num : ('d, 'k) session -> int * int -> unit
      = "sunml_arkode_ark_set_table_num"

  let set_erk_table_num s v =
    c_set_table_num s (-1, ButcherTable.int_of_erk_table v)

  let set_dirk_table_num s v =
    c_set_table_num s (ButcherTable.int_of_dirk_table v, -1)

  let set_ark_table_num s v =
    c_set_table_num s (ButcherTable.ints_of_ark_table v)

  external c_set_adaptivity_method
      : ('d, 'k) session -> 'd adaptivity_method -> unit
      = "sunml_arkode_ark_set_adaptivity_method"

  let set_adaptivity_method s am =
    (match am with
     | AdaptivityFn fn -> s.adaptfn <- fn
     | _ -> s.adaptfn <- dummy_adaptfn);
    c_set_adaptivity_method s am

  external c_set_stability_fn : ('d, 'k) session -> bool -> unit
      = "sunml_arkode_ark_set_stability_fn"

  let set_stability_fn s f =
    s.stabfn <- f;
    c_set_stability_fn s true

  let clear_stability_fn s =
    s.stabfn <- dummy_stabfn;
    c_set_stability_fn s false

  type 'd stage_predict_fn  = 'd Arkode_impl.Global.stage_predict_fn

  external set_predictor_method : ('d, 'k) session -> predictor_method -> unit
      = "sunml_arkode_ark_set_predictor_method"

  external c_set_stage_predict_fn
      : ('d, 'k) session -> bool -> unit
      = "sunml_arkode_ark_set_stage_predict_fn"

  let set_stage_predict_fn s fn =
    s.stagepredictfn <- fn;
    c_set_stage_predict_fn s true

  let clear_stage_predict_fn s =
    c_set_stage_predict_fn s false;
    s.stagepredictfn <- dummy_stagepredictfn

  external set_defaults           : ('a, 'k) session -> unit
      = "sunml_arkode_ark_set_defaults"
  external set_interpolant_type : ('d, 'k) session -> interpolant_type -> unit
      = "sunml_arkode_ark_set_interpolant_type"
  external set_interpolant_degree : ('d, 'k) session -> int -> unit
      = "sunml_arkode_ark_set_interpolant_degree"
  external set_max_num_steps      : ('a, 'k) session -> int -> unit
      = "sunml_arkode_ark_set_max_num_steps"
  external set_max_hnil_warns     : ('a, 'k) session -> int -> unit
      = "sunml_arkode_ark_set_max_hnil_warns"
  external set_init_step          : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_init_step"
  external c_set_fixed_step         : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_fixed_step"
  let set_fixed_step s ohf =
    c_set_fixed_step s (match ohf with None -> 0.0 | Some v -> v)
  external set_max_num_constr_fails : ('a, 'k) session -> int -> unit
      = "sunml_arkode_ark_set_max_num_constr_fails"
  external set_min_step           : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_min_step"
  external set_max_step           : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_max_step"
  external set_stop_time          : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_stop_time"
  external set_optimal_params     : ('a, 'k) session -> unit
      = "sunml_arkode_ark_set_optimal_params"
  external set_max_err_test_fails : ('a, 'k) session -> int -> unit
      = "sunml_arkode_ark_set_max_err_test_fails"
  external set_max_nonlin_iters   : ('a, 'k) session -> int -> unit
      = "sunml_arkode_ark_set_max_nonlin_iters"
  external set_max_conv_fails     : ('a, 'k) session -> int -> unit
      = "sunml_arkode_ark_set_max_conv_fails"
  external set_nonlin_conv_coef   : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_nonlin_conv_coef"
  external set_constraints      : ('a, 'k) session -> ('a, 'k) Nvector.t -> unit
      = "sunml_arkode_ark_set_constraints"
  external set_nonlin_crdown      : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_nonlin_crdown"
  external set_nonlin_rdiv        : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_nonlin_rdiv"
  external set_delta_gamma_max    : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_delta_gamma_max"
  external set_lsetup_frequency   : ('a, 'k) session -> int -> unit
      = "sunml_arkode_ark_set_lsetup_frequency"
  external set_cfl_fraction       : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_cfl_fraction"
  external set_error_bias         : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_error_bias"
  external set_fixed_step_bounds  : ('a, 'k) session -> float -> float -> unit
      = "sunml_arkode_ark_set_fixed_step_bounds"
  external set_max_cfail_growth   : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_max_cfail_growth"
  external set_max_efail_growth   : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_max_efail_growth"
  external set_max_first_growth   : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_max_first_growth"
  external set_max_growth         : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_max_growth"
  external set_min_reduction      : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_min_reduction"
  external set_safety_factor      : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_safety_factor"
  external set_small_num_efails   : ('a, 'k) session -> float -> unit
      = "sunml_arkode_ark_set_small_num_efails"

  external c_set_postprocess_step_fn : ('a, 'k) session -> bool -> unit
      = "sunml_arkode_ark_set_postprocess_step_fn"

  let set_postprocess_step_fn s fn =
    s.poststepfn <- fn;
    c_set_postprocess_step_fn s true

  let clear_postprocess_step_fn s =
    c_set_postprocess_step_fn s false;
    s.poststepfn <- dummy_poststepfn

  external c_set_root_direction   : ('a, 'k) session -> RootDirs.t -> unit
      = "sunml_arkode_ark_set_root_direction"

  let set_root_direction s rda =
    c_set_root_direction s (RootDirs.copy (get_num_roots s) rda)

  let set_all_root_directions s rd =
    c_set_root_direction s (RootDirs.make (get_num_roots s) rd)

  external set_no_inactive_root_warn      : ('a, 'k) session -> unit
      = "sunml_arkode_ark_set_no_inactive_root_warn"

  external get_current_butcher_tables
      : ('d, 'k) session -> ButcherTable.t option * ButcherTable.t option
      = "sunml_arkode_ark_get_current_butcher_tables"

  external get_tol_scale_factor           : ('a, 'k) session -> float
      = "sunml_arkode_ark_get_tol_scale_factor"

  external c_get_err_weights
      : ('a, 'k) session -> ('a, 'k) nvector -> unit
      = "sunml_arkode_ark_get_err_weights"

  let get_err_weights s ew =
    if Sundials_configuration.safe then s.checkvec ew;
    c_get_err_weights s ew

  external c_get_res_weights
      : ('a, 'k) session -> ('a, 'k) nvector -> unit
      = "sunml_arkode_ark_get_res_weights"

  let get_res_weights s rw =
    if Sundials_configuration.safe then s.checkvec rw;
    c_get_res_weights s rw

  external c_get_est_local_errors
      : ('a, 'k) session -> ('a, 'k) nvector -> unit
      = "sunml_arkode_ark_get_est_local_errors"

  let get_est_local_errors s ew =
    if Sundials_configuration.safe then s.checkvec ew;
    c_get_est_local_errors s ew

  external get_num_nonlin_solv_iters      : ('a, 'k) session -> int
      = "sunml_arkode_ark_get_num_nonlin_solv_iters"

  external get_num_nonlin_solv_conv_fails : ('a, 'k) session -> int
      = "sunml_arkode_ark_get_num_nonlin_solv_conv_fails"

  external get_nonlin_solv_stats          : ('a, 'k) session -> int * int
      = "sunml_arkode_ark_get_nonlin_solv_stats"

  external get_num_g_evals                : ('a, 'k) session -> int
      = "sunml_arkode_ark_get_num_g_evals"

  external get_num_constr_fails           : ('a, 'k) session -> int
      = "sunml_arkode_ark_get_num_constr_fails"

  external write_parameters : ('d, 'k) session -> Logfile.t -> unit
      = "sunml_arkode_ark_write_parameters"

  external write_butcher : ('d, 'k) session -> Logfile.t -> unit
      = "sunml_arkode_ark_write_butcher"

end (* }}} *)

module ERKStep = struct (* {{{ *)
  include Arkode_impl
  include Common

  type ('d, 'k) session = ('d, 'k, erkstep) Arkode_impl.session

  external c_root_init : ('a, 'k) session -> int -> unit
      = "sunml_arkode_erk_root_init"

  let root_init session (nroots, rootsfn) =
    c_root_init session nroots;
    session.rootsfn <- rootsfn

  external sv_tolerances
      : ('a, 'k) session -> float -> ('a, 'k) nvector -> unit
      = "sunml_arkode_erk_sv_tolerances"
  external ss_tolerances
      : ('a, 'k) session -> float -> float -> unit
      = "sunml_arkode_erk_ss_tolerances"
  external wf_tolerances 
      : ('a, 'k) session -> unit
      = "sunml_arkode_erk_wf_tolerances"

  let set_tolerances s tol =
    match tol with
    | SStolerances (rel, abs) -> (s.errw <- dummy_errw; ss_tolerances s rel abs)
    | SVtolerances (rel, abs) -> (if Sundials_configuration.safe then s.checkvec abs;
                                  s.errw <- dummy_errw; sv_tolerances s rel abs)
    | WFtolerances ferrw -> (s.errw <- ferrw; wf_tolerances s)

  external c_set_order : ('a, 'k) session -> int -> unit
    = "sunml_arkode_erk_set_order"

  external session_finalize : ('a, 'k) session -> unit
      = "sunml_arkode_erk_session_finalize"

  external c_init :
    ('a, 'k) session Weak.t
    -> ('a, 'k) nvector (* y_0 *)
    -> float            (* t_0 *)
    -> (erkstep arkode_mem * c_weak_ref)
    = "sunml_arkode_erk_init"

  let init tol ?order f ?(roots=no_roots) t0 y0 =
    let (nroots, roots) = roots in
    let checkvec = Nvector.check y0 in
    if Sundials_configuration.safe && nroots < 0 then
      raise (Invalid_argument "number of root functions is negative");
    let weakref = Weak.create 1 in
    let arkode_mem, backref = c_init weakref y0 t0 in
    (* arkode_mem and backref have to be immediately captured in a session and
       associated with the finalizer before we do anything else.  *)
    let session = {
            arkode       = arkode_mem;
            backref      = backref;
            nroots       = nroots;
            checkvec     = checkvec;
            uses_resv    = false;

            exn_temp     = None;

            problem      = ExplicitOnly;
            rhsfn1       = f;
            rhsfn2       = dummy_rhsfn2;

            rootsfn      = roots;
            errh         = dummy_errh;
            errw         = dummy_errw;
            resw         = dummy_resw;

            adaptfn      = dummy_adaptfn;
            stabfn       = dummy_stabfn;
            resizefn     = dummy_resizefn;
            poststepfn   = dummy_poststepfn;
            stagepredictfn = dummy_stagepredictfn;
            preinnerfn   = dummy_preinnerfn;
            postinnerfn  = dummy_postinnerfn;
            preinnerarray = empty_preinnerarray ();

            linsolver      = None;
            ls_solver      = LSI.NoHLS;
            ls_callbacks   = NoCallbacks;
            ls_precfns     = NoPrecFns;

            mass_solver    = LSI.NoHLS;
            mass_callbacks = NoMassCallbacks;
            mass_precfns   = NoMassPrecFns;

            nls_solver     = None;

            inner_session  = None;
          } in
    Gc.finalise session_finalize session;
    Weak.set weakref 0 (Some session);
    (* Now the session is safe to use.  If any of the following fails and raises
       an exception, the GC will take care of freeing arkode_mem and backref.  *)
    if nroots > 0 then
      c_root_init session nroots;
    set_tolerances session tol;
    (match order with Some o -> c_set_order session o | None -> ());
    session

  let get_num_roots { nroots } = nroots

  external reset : ('d, 'k) session -> float -> ('d, 'k) Nvector.t
      = "sunml_arkode_erk_reset"

  external c_reinit
      : ('a, 'k) session -> float -> ('a, 'k) nvector -> unit
      = "sunml_arkode_erk_reinit"

  let reinit session ?order ?roots t0 y0 =
    if Sundials_configuration.safe then session.checkvec y0;
    c_reinit session t0 y0;
    (match order with Some o -> c_set_order session o | None -> ());
    (match roots with Some roots -> root_init session roots| None -> ())

  external c_resize
      : ('a, 'k) session -> bool -> float -> float -> ('a, 'k) nvector -> unit
      = "sunml_arkode_erk_resize"

  let resize session ?resize_nvec tol hscale ynew t0 =
    session.checkvec <- Nvector.check ynew;
    (match resize_nvec with None -> () | Some f -> session.resizefn <- f);
    c_resize session (resize_nvec <> None) hscale t0 ynew;
    session.resizefn <- dummy_resizefn;
    set_tolerances session tol

  external get_root_info  : ('a, 'k) session -> Roots.t -> unit
      = "sunml_arkode_erk_get_root_info"

  external c_solve_normal : ('a, 'k) session -> float -> ('a, 'k) nvector
                                -> float * solver_result
      = "sunml_arkode_erk_solve_normal"

  let solve_normal s t y =
    if Sundials_configuration.safe then s.checkvec y;
    c_solve_normal s t y

  external c_solve_one_step : ('a, 'k) session -> float -> ('a, 'k) nvector
                                -> float * solver_result
      = "sunml_arkode_erk_solve_one_step"

  let solve_one_step s t y =
    if Sundials_configuration.safe then s.checkvec y;
    c_solve_one_step s t y

  external c_get_dky
      : ('a, 'k) session -> float -> int -> ('a, 'k) nvector -> unit
      = "sunml_arkode_erk_get_dky"

  let get_dky s y =
    if Sundials_configuration.safe then s.checkvec y;
    fun t k -> c_get_dky s t k y

  (* Synchronized with arkode_timestepper_stats_index in arkode_ml.h *)
  type timestepper_stats = {
      exp_steps           : int;
      acc_steps           : int;
      step_attempts       : int;
      num_nf_evals        : int;
      num_err_test_fails  : int;
    }

  external get_timestepper_stats : ('a, 'k) session -> timestepper_stats
      = "sunml_arkode_erk_get_timestepper_stats"

  external get_step_stats : ('a, 'k) session -> step_stats
      = "sunml_arkode_erk_get_step_stats"

  external get_work_space         : ('a, 'k) session -> int * int
      = "sunml_arkode_erk_get_work_space"

  external get_num_steps          : ('a, 'k) session -> int
      = "sunml_arkode_erk_get_num_steps"

  external get_num_exp_steps      : ('d, 'k) session -> int
      = "sunml_arkode_erk_get_num_exp_steps"

  external get_num_acc_steps      : ('d, 'k) session -> int
      = "sunml_arkode_erk_get_num_acc_steps"

  external get_num_step_attempts  : ('d, 'k) session -> int
      = "sunml_arkode_erk_get_num_step_attempts"

  external get_num_rhs_evals      : ('a, 'k) session -> int
      = "sunml_arkode_erk_get_num_rhs_evals"

  external get_num_err_test_fails : ('a, 'k) session -> int
      = "sunml_arkode_erk_get_num_err_test_fails"

  external get_actual_init_step   : ('a, 'k) session -> float
      = "sunml_arkode_erk_get_actual_init_step"

  external get_last_step          : ('a, 'k) session -> float
      = "sunml_arkode_erk_get_last_step"

  external get_current_step       : ('a, 'k) session -> float
      = "sunml_arkode_erk_get_current_step"

  external get_current_time       : ('a, 'k) session -> float
      = "sunml_arkode_erk_get_current_time"

  let print_timestepper_stats s oc =
    let stats = get_timestepper_stats s
    in
    Printf.fprintf oc "exp_steps = %d\n"           stats.exp_steps;
    Printf.fprintf oc "acc_steps = %d\n"           stats.acc_steps;
    Printf.fprintf oc "step_attempts = %d\n"       stats.step_attempts;
    Printf.fprintf oc "num_nf_evals = %d\n"        stats.num_nf_evals;
    Printf.fprintf oc "num_err_test_fails = %d\n"  stats.num_err_test_fails

  let print_step_stats s oc =
    print_step_stats oc (get_step_stats s)

  external set_diagnostics : ('a, 'k) session -> Logfile.t -> unit
      = "sunml_arkode_erk_set_diagnostics"

  external clear_diagnostics : ('a, 'k) session -> unit
      = "sunml_arkode_erk_clear_diagnostics"

  external set_interpolant_type : ('d, 'k) session -> interpolant_type -> unit
      = "sunml_arkode_erk_set_interpolant_type"

  external set_interpolant_degree : ('d, 'k) session -> int -> unit
      = "sunml_arkode_erk_set_interpolant_degree"

  external set_error_file : ('a, 'k) session -> Logfile.t -> unit
      = "sunml_arkode_erk_set_error_file"

  external c_set_err_handler_fn  : ('a, 'k) session -> unit
      = "sunml_arkode_erk_set_err_handler_fn"

  let set_err_handler_fn s ferrh =
    s.errh <- ferrh;
    c_set_err_handler_fn s

  external clear_err_handler_fn  : ('a, 'k) session -> unit
      = "sunml_arkode_erk_clear_err_handler_fn"

  let clear_err_handler_fn s =
    s.errh <- dummy_errh;
    clear_err_handler_fn s

  external c_set_table
    : ('d, 'k) session -> ButcherTable.t option -> unit
    = "sunml_arkode_erk_set_table"

  let set_table s table =
    c_set_table s (Some table)

  external c_set_table_num : ('d, 'k) session -> int -> unit
      = "sunml_arkode_erk_set_table_num"

  let set_table_num s v =
    c_set_table_num s (ButcherTable.int_of_erk_table v)

  external c_set_adaptivity_method
      : ('d, 'k) session -> 'd adaptivity_method -> unit
      = "sunml_arkode_erk_set_adaptivity_method"

  let set_adaptivity_method s am =
    (match am with
     | AdaptivityFn fn -> s.adaptfn <- fn
     | _ -> s.adaptfn <- dummy_adaptfn);
    c_set_adaptivity_method s am

  external c_set_stability_fn : ('d, 'k) session -> bool -> unit
      = "sunml_arkode_erk_set_stability_fn"

  let set_stability_fn s f =
    s.stabfn <- f;
    c_set_stability_fn s true

  let clear_stability_fn s =
    s.stabfn <- dummy_stabfn;
    c_set_stability_fn s false

  external set_defaults           : ('a, 'k) session -> unit
      = "sunml_arkode_erk_set_defaults"
  external c_set_fixed_step       : ('a, 'k) session -> float -> unit
      = "sunml_arkode_erk_set_fixed_step"
  let set_fixed_step s ohf =
    c_set_fixed_step s (match ohf with None -> 0.0 | Some v -> v)
  external set_max_num_constr_fails : ('a, 'k) session -> int -> unit
      = "sunml_arkode_erk_set_max_num_constr_fails"
  external set_init_step          : ('a, 'k) session -> float -> unit
      = "sunml_arkode_erk_set_init_step"
  external set_max_hnil_warns     : ('a, 'k) session -> int -> unit
      = "sunml_arkode_erk_set_max_hnil_warns"
  external set_max_num_steps      : ('a, 'k) session -> int -> unit
      = "sunml_arkode_erk_set_max_num_steps"
  external set_max_step           : ('a, 'k) session -> float -> unit
      = "sunml_arkode_erk_set_max_step"
  external set_min_step           : ('a, 'k) session -> float -> unit
      = "sunml_arkode_erk_set_min_step"
  external set_stop_time          : ('a, 'k) session -> float -> unit
      = "sunml_arkode_erk_set_stop_time"
  external set_max_err_test_fails : ('a, 'k) session -> int -> unit
      = "sunml_arkode_erk_set_max_err_test_fails"
  external set_cfl_fraction       : ('a, 'k) session -> float -> unit
      = "sunml_arkode_erk_set_cfl_fraction"
  external set_error_bias         : ('a, 'k) session -> float -> unit
      = "sunml_arkode_erk_set_error_bias"
  external set_fixed_step_bounds  : ('a, 'k) session -> float -> float -> unit
      = "sunml_arkode_erk_set_fixed_step_bounds"
  external set_max_efail_growth   : ('a, 'k) session -> float -> unit
      = "sunml_arkode_erk_set_max_efail_growth"
  external set_max_first_growth   : ('a, 'k) session -> float -> unit
      = "sunml_arkode_erk_set_max_first_growth"
  external set_max_growth         : ('a, 'k) session -> float -> unit
      = "sunml_arkode_erk_set_max_growth"
  external set_min_reduction      : ('a, 'k) session -> float -> unit
      = "sunml_arkode_erk_set_min_reduction"
  external set_safety_factor      : ('a, 'k) session -> float -> unit
      = "sunml_arkode_erk_set_safety_factor"
  external set_small_num_efails   : ('a, 'k) session -> float -> unit
      = "sunml_arkode_erk_set_small_num_efails"
  external set_constraints      : ('a, 'k) session -> ('a, 'k) Nvector.t -> unit
      = "sunml_arkode_erk_set_constraints"

  external c_set_postprocess_step_fn : ('a, 'k) session -> bool -> unit
      = "sunml_arkode_erk_set_postprocess_step_fn"

  let set_postprocess_step_fn s fn =
    s.poststepfn <- fn;
    c_set_postprocess_step_fn s true

  let clear_postprocess_step_fn s =
    s.poststepfn <- dummy_poststepfn;
    c_set_postprocess_step_fn s false

  external c_set_root_direction   : ('a, 'k) session -> RootDirs.t -> unit
      = "sunml_arkode_erk_set_root_direction"

  let set_root_direction s rda =
    c_set_root_direction s (RootDirs.copy (get_num_roots s) rda)

  let set_all_root_directions s rd =
    c_set_root_direction s (RootDirs.make (get_num_roots s) rd)

  external set_no_inactive_root_warn      : ('a, 'k) session -> unit
      = "sunml_arkode_erk_set_no_inactive_root_warn"

  external get_current_butcher_table
      : ('d, 'k) session -> ButcherTable.t
      = "sunml_arkode_erk_get_current_butcher_table"

  external get_tol_scale_factor           : ('a, 'k) session -> float
      = "sunml_arkode_erk_get_tol_scale_factor"

  external c_get_err_weights
      : ('a, 'k) session -> ('a, 'k) nvector -> unit
      = "sunml_arkode_erk_get_err_weights"

  let get_err_weights s ew =
    if Sundials_configuration.safe then s.checkvec ew;
    c_get_err_weights s ew

  external c_get_est_local_errors
      : ('a, 'k) session -> ('a, 'k) nvector -> unit
      = "sunml_arkode_erk_get_est_local_errors"

  let get_est_local_errors s ew =
    if Sundials_configuration.safe then s.checkvec ew;
    c_get_est_local_errors s ew

  external get_num_g_evals                : ('a, 'k) session -> int
      = "sunml_arkode_erk_get_num_g_evals"

  external get_num_constr_fails           : ('a, 'k) session -> int
      = "sunml_arkode_erk_get_num_constr_fails"

  external write_parameters : ('d, 'k) session -> Logfile.t -> unit
      = "sunml_arkode_erk_write_parameters"

  external write_butcher : ('d, 'k) session -> Logfile.t -> unit
      = "sunml_arkode_erk_write_butcher"

end (* }}} *)

module MRIStep = struct (* {{{ *)
  include Arkode_impl
  include Common

  type ('d, 'k) session = ('d, 'k, mristep) Arkode_impl.session

  type 'k serial_session = (Nvector_serial.data, 'k) session
                           constraint 'k = [>Nvector_serial.kind]

  type ('d, 'k) linear_solver = ('d, 'k, mristep) lin_solver

  type 'k serial_linear_solver =
    (Nvector_serial.data, 'k) linear_solver
    constraint 'k = [>Nvector_serial.kind]

  external c_root_init : ('a, 'k) session -> int -> unit
      = "sunml_arkode_mri_root_init"

  let root_init session (nroots, rootsfn) =
    c_root_init session nroots;
    session.rootsfn <- rootsfn

  module Dls = struct (* {{{ *)
    include Dls
    include LinearSolver.Direct

    let check_dqjac (type k m nd nk) jac (mat : (k,m,nd,nk) Matrix.t) =
      let open Matrix in
      match get_id mat with
      | Dense | Band -> ()
      | _ -> if jac = None then invalid_arg "A Jacobian function is required"

    let set_ls_callbacks (type m) (type tag)
          ?jac ?(linsys : m linsys_fn option)
          (solver_data : (m, 'nd, 'nk, tag) LSI.solver_data)
          (mat : ('mk, m, 'nd, 'nk) Matrix.t) session =
      let cb = { jacfn = (match jac with None -> no_callback | Some f -> f);
                 jmat  = (None : m option) } in
      let ls = match linsys with None -> DirectTypes.no_linsysfn | Some f -> f in
      begin match solver_data with
      | LSI.Dense ->
          session.ls_callbacks <- DlsDenseCallback (cb, ls)
      | LSI.LapackDense ->
          session.ls_callbacks <- DlsDenseCallback (cb, ls)
      | LSI.Band ->
          session.ls_callbacks <- DlsBandCallback (cb, ls)
      | LSI.LapackBand ->
          session.ls_callbacks <- DlsBandCallback (cb, ls)
      | LSI.Klu _ ->
          if jac = None then invalid_arg "Klu requires Jacobian function";
          session.ls_callbacks <- SlsKluCallback (cb, ls)
      | LSI.Superlumt _ ->
          if jac = None then invalid_arg "Superlumt requires Jacobian function";
          session.ls_callbacks <- SlsSuperlumtCallback (cb, ls)
      | LSI.Custom _ ->
          check_dqjac jac mat;
          session.ls_callbacks <- DirectCustomCallback (cb, ls)
      | _ -> assert false
      end;
      session.ls_precfns <- NoPrecFns

    external c_set_linear_solver
      : ('d, 'k) session
        -> ('m, 'd, 'k) LSI.cptr
        -> ('mk, 'm, 'd, 'k) Matrix.t option
        -> bool
        -> bool
        -> unit
      = "sunml_arkode_mri_set_linear_solver"

    let assert_matrix = function
      | Some m -> m
      | None -> failwith "a direct linear solver is required"

    let solver ?jac ?linsys ls session nv =
      if Sundials_impl.Version.lt540
        then raise Config.NotImplementedBySundialsVersion;
      let LSI.LS ({ rawptr; solver; matrix } as hls) = ls in
      let matrix = assert_matrix matrix in
      set_ls_callbacks ?jac ?linsys solver matrix session;
      c_set_linear_solver session rawptr (Some matrix) (jac <> None)
                                                       (linsys <> None);
      LSI.attach ls;
      session.ls_solver <- LSI.HLS hls

    external get_work_space : 'k serial_session -> int * int
      = "sunml_arkode_mri_get_lin_work_space"

    external get_num_jac_evals : 'k serial_session -> int
      = "sunml_arkode_mri_get_num_jac_evals"

    external get_num_lin_rhs_evals : 'k serial_session -> int
      = "sunml_arkode_mri_get_num_lin_rhs_evals"

  end (* }}} *)

  module Spils = struct (* {{{ *)
    include Spils
    include LinearSolver.Iterative

    external c_set_jac_times : ('a, 'k) session -> bool -> bool -> unit
      = "sunml_arkode_mri_set_jac_times"

    external c_set_jac_times_rhsfn : ('a, 'k) session -> bool -> unit
      = "sunml_arkode_mri_set_jac_times_rhsfn"

    external c_set_preconditioner
      : ('a, 'k) session -> bool -> unit
      = "sunml_arkode_mri_set_preconditioner"

    let init_preconditioner solve setup session nv =
      c_set_preconditioner session (setup <> None);
      session.ls_precfns <- PrecFns { prec_solve_fn = solve;
                                      prec_setup_fn = setup }

    let prec_none = LSI.Iterative.(PrecNone,
                      fun session nv -> session.ls_precfns <- NoPrecFns)

    let prec_left ?setup solve  = LSI.Iterative.(PrecLeft,
                                              init_preconditioner solve setup)

    let prec_right ?setup solve = LSI.Iterative.(PrecRight,
                                              init_preconditioner solve setup)

    let prec_both ?setup solve  = LSI.Iterative.(PrecBoth,
                                              init_preconditioner solve setup)

    external c_set_linear_solver
      : ('d, 'k) session
        -> ('m, 'd, 'k) LSI.cptr
        -> ('mk, 'm, 'd, 'k) Matrix.t option
        -> bool
        -> bool
        -> unit
      = "sunml_arkode_ark_set_linear_solver"

    let solver (type s)
          (LSI.(LS ({ rawptr; solver; compat } as hls) as ls))
          ?jac_times_vec ?jac_times_rhs (prec_type, set_prec) session nv =
      if Sundials_impl.Version.lt540
        then raise Config.NotImplementedBySundialsVersion;
      let jac_times_setup, jac_times_vec =
        match jac_times_vec with
        | None -> None, None
        | Some _ when jac_times_rhs <> None ->
            invalid_arg "cannot pass both jac_times_vec and jac_times_rhs"
        | Some (ojts, jtv) -> ojts, Some jtv
      in
      c_set_linear_solver session rawptr None false false;
      LSI.attach ls;
      session.ls_solver <- LSI.HLS hls;
      LSI.(impl_set_prec_type rawptr solver prec_type false);
      set_prec session nv;
      match jac_times_rhs with
      | Some jtrhsfn -> begin
          session.ls_callbacks <- SpilsCallback2 jtrhsfn;
          c_set_jac_times_rhsfn session true
        end
      | None -> begin
          session.ls_callbacks <-
            SpilsCallback1 (jac_times_vec, jac_times_setup);
          if jac_times_setup <> None || jac_times_vec <> None then
            c_set_jac_times session (jac_times_setup <> None)
                                    (jac_times_vec <> None)
        end

    let set_jac_times s ?jac_times_setup f =
      if Sundials_impl.Version.lt540
        then raise Config.NotImplementedBySundialsVersion;
      match s.ls_callbacks with
      | SpilsCallback1 _ ->
          c_set_jac_times s (jac_times_setup <> None) true;
          s.ls_callbacks <- SpilsCallback1 (Some f, jac_times_setup)
      | _ -> raise LinearSolver.InvalidLinearSolver

    let clear_jac_times s =
      if Sundials_impl.Version.lt540
        then raise Config.NotImplementedBySundialsVersion;
      match s.ls_callbacks with
      | SpilsCallback1 _ ->
          c_set_jac_times s false false;
          s.ls_callbacks <- SpilsCallback1 (None, None)
      | _ -> raise LinearSolver.InvalidLinearSolver

    let set_preconditioner s ?setup solve =
      if Sundials_impl.Version.lt540
        then raise Config.NotImplementedBySundialsVersion;
      match s.ls_callbacks with
      | SpilsCallback1 _ | SpilsCallback2 _ ->
          c_set_preconditioner s (setup <> None);
          s.ls_precfns <- PrecFns { prec_setup_fn = setup;
                                    prec_solve_fn = solve }
      | _ -> raise LinearSolver.InvalidLinearSolver

    external set_jac_eval_frequency : ('a, 'k) session -> int -> unit
        = "sunml_arkode_mri_set_jac_eval_frequency"

    external set_linear_solution_scaling : ('d, 'k) session -> bool -> unit
      = "sunml_arkode_mri_set_linear_solution_scaling"

    external set_eps_lin            : ('a, 'k) session -> float -> unit
      = "sunml_arkode_mri_set_eps_lin"

    external set_ls_norm_factor : ('d, 'k) session -> float -> unit
      = "sunml_arkode_mri_set_ls_norm_factor"

    external get_num_lin_iters      : ('a, 'k) session -> int
      = "sunml_arkode_mri_get_num_lin_iters"

    external get_num_lin_conv_fails : ('a, 'k) session -> int
      = "sunml_arkode_mri_get_num_lin_conv_fails"

    external get_work_space         : ('a, 'k) session -> int * int
      = "sunml_arkode_mri_get_lin_work_space"

    external get_num_prec_evals     : ('a, 'k) session -> int
      = "sunml_arkode_mri_get_num_prec_evals"

    external get_num_prec_solves    : ('a, 'k) session -> int
      = "sunml_arkode_mri_get_num_prec_solves"

    let get_num_prec_solves s =
      ls_check_spils s;
      get_num_prec_solves s

    external get_num_jtsetup_evals   : ('a, 'k) session -> int
      = "sunml_arkode_mri_get_num_jtsetup_evals"

    let get_num_jtsetup_evals s =
      ls_check_spils s;
      get_num_jtsetup_evals s

    external get_num_jtimes_evals   : ('a, 'k) session -> int
      = "sunml_arkode_mri_get_num_jtimes_evals"

    let get_num_jtimes_evals s =
      ls_check_spils s;
      get_num_jtimes_evals s

    external get_num_lin_rhs_evals  : ('a, 'k) session -> int
      = "sunml_arkode_mri_get_num_lin_rhs_evals"

    let get_num_lin_rhs_evals s =
      ls_check_spils s;
      get_num_lin_rhs_evals s

  end (* }}} *)

  external session_finalize : ('a, 'k) session -> unit
      = "sunml_arkode_mri_session_finalize"

  external c_set_fixed_step : ('a, 'k) session -> float -> unit
      = "sunml_arkode_mri_set_fixed_step"

  external sv_tolerances
      : ('a, 'k) session -> float -> ('a, 'k) nvector -> unit
      = "sunml_arkode_mri_sv_tolerances"
  external ss_tolerances
      : ('a, 'k) session -> float -> float -> unit
      = "sunml_arkode_mri_ss_tolerances"
  external wf_tolerances 
      : ('a, 'k) session -> unit
      = "sunml_arkode_mri_wf_tolerances"

  let set_tolerances s tol =
    if Sundials_impl.Version.lt540 && tol <> default_tolerances
      then raise Config.NotImplementedBySundialsVersion;
    match tol with
    | SStolerances (rel, abs) -> (s.errw <- dummy_errw; ss_tolerances s rel abs)
    | SVtolerances (rel, abs) -> (if Sundials_configuration.safe then s.checkvec abs;
                                  s.errw <- dummy_errw; sv_tolerances s rel abs)
    | WFtolerances ferrw -> (s.errw <- ferrw; wf_tolerances s)

  external c_init :
    ('a, 'k) session Weak.t
    -> ('a, 'k) ARKStep.session
    -> ('a, 'k) nvector (* y_0 *)
    -> float            (* t_0 *)
    -> (mristep arkode_mem * c_weak_ref)
    = "sunml_arkode_mri_init"

  (* 5.4.0 <= Sundials *)
  external c_set_nonlinear_solver
      : ('d, 'k) session
        -> ('d, 'k, ('d, 'k) session, [`Nvec]) NLSI.cptr
        -> unit
      = "sunml_arkode_mri_set_nonlinear_solver"

  external set_linear : ('a, 'k) session -> bool -> unit
    = "sunml_arkode_mri_set_linear"

  let init fasts tol ?nlsolver ?lsolver ?linearity slow
           ~slowstep ?(roots=no_roots) t0 y0 =
    if Sundials_impl.Version.lt500
      then raise Config.NotImplementedBySundialsVersion;
    let (nroots, roots) = roots in
    let checkvec = Nvector.check y0 in
    if Sundials_configuration.safe && nroots < 0
      then raise Config.NotImplementedBySundialsVersion;
    if Sundials_impl.Version.lt540
        && (nlsolver <> None || lsolver <> None || linearity <> None)
      then invalid_arg "functions is negative";
    let weakref = Weak.create 1 in
    let arkode_mem, backref = c_init weakref fasts y0 t0 in
    (* arkode_mem and backref have to be immediately captured in a session and
       associated with the finalizer before we do anything else.  *)
    let session = {
            arkode       = arkode_mem;
            backref      = backref;
            nroots       = nroots;
            checkvec     = checkvec;
            uses_resv    = false;

            exn_temp     = None;

            problem      = ExplicitOnly; (* ignored for MRIStep *)
            rhsfn1       = slow;
            rhsfn2       = dummy_rhsfn2;

            rootsfn      = roots;
            errh         = dummy_errh;
            errw         = dummy_errw;
            resw         = dummy_resw;

            adaptfn      = dummy_adaptfn;
            stabfn       = dummy_stabfn;
            resizefn     = dummy_resizefn;
            poststepfn   = dummy_poststepfn;
            stagepredictfn = dummy_stagepredictfn;
            preinnerfn   = dummy_preinnerfn;
            postinnerfn  = dummy_postinnerfn;
            preinnerarray = empty_preinnerarray ();

            linsolver      = None;
            ls_solver      = LSI.NoHLS;
            ls_callbacks   = NoCallbacks;
            ls_precfns     = NoPrecFns;

            mass_solver    = LSI.NoHLS;
            mass_callbacks = NoMassCallbacks;
            mass_precfns   = NoMassPrecFns;

            nls_solver     = None;

            inner_session  = Some fasts;
          } in
    Gc.finalise session_finalize session;
    Weak.set weakref 0 (Some session);
    (* Now the session is safe to use.  If any of the following fails and raises
       an exception, the GC will take care of freeing arkode_mem and backref.  *)
    if nroots > 0 then c_root_init session nroots;
    set_tolerances session tol;
    (match nlsolver with
     | Some ({ NLSI.rawptr = nlcptr } as nls) ->
         NLSI.attach nls;
         session.nls_solver <- Some nls;
         c_set_nonlinear_solver session nlcptr
     | _ -> ());
    (match lsolver with
     | None -> ()
     | Some ls -> session.linsolver <- Some ls;
                  ls session y0);
    (match linearity with
     | Some (Linear timedepend) -> set_linear session timedepend
     | _ -> ());
    c_set_fixed_step session slowstep;
    session

  let get_num_roots { nroots } = nroots

  external reset : ('d, 'k) session -> float -> ('d, 'k) Nvector.t
      = "sunml_arkode_mri_reset"

  external c_reinit
      : ('a, 'k) session -> float -> ('a, 'k) nvector -> unit
      = "sunml_arkode_mri_reinit"

  let reinit session ?nlsolver ?lsolver ?roots t0 y0 =
    if Sundials_configuration.safe then session.checkvec y0;
    c_reinit session t0 y0;
    (match roots with Some roots -> root_init session roots| None -> ())

  external c_resize
      : ('a, 'k) session -> bool -> float -> ('a, 'k) nvector -> unit
      = "sunml_arkode_mri_resize"

  let resize session ?resize_nvec ynew t0 =
    session.checkvec <- Nvector.check ynew;
    (match resize_nvec with None -> () | Some f -> session.resizefn <- f);
    c_resize session (resize_nvec <> None) t0 ynew;
    session.resizefn <- dummy_resizefn

  external get_root_info  : ('a, 'k) session -> Roots.t -> unit
      = "sunml_arkode_mri_get_root_info"

  external c_solve_normal : ('a, 'k) session -> float -> ('a, 'k) nvector
                                -> float * solver_result
      = "sunml_arkode_mri_solve_normal"

  let solve_normal s t y =
    if Sundials_configuration.safe then s.checkvec y;
    c_solve_normal s t y

  external c_solve_one_step : ('a, 'k) session -> float -> ('a, 'k) nvector
                                -> float * solver_result
      = "sunml_arkode_mri_solve_one_step"

  let solve_one_step s t y =
    if Sundials_configuration.safe then s.checkvec y;
    c_solve_one_step s t y

  external c_get_dky
      : ('a, 'k) session -> float -> int -> ('a, 'k) nvector -> unit
      = "sunml_arkode_mri_get_dky"

  let get_dky s y =
    if Sundials_configuration.safe then s.checkvec y;
    fun t k -> c_get_dky s t k y

  external get_work_space         : ('a, 'k) session -> int * int
      = "sunml_arkode_mri_get_work_space"

  external get_num_steps          : ('a, 'k) session -> int
      = "sunml_arkode_mri_get_num_steps"

  external get_num_rhs_evals      : ('a, 'k) session -> int
      = "sunml_arkode_mri_get_num_rhs_evals"

  external get_last_step          : ('a, 'k) session -> float
      = "sunml_arkode_mri_get_last_step"

  external get_current_time       : ('a, 'k) session -> float
      = "sunml_arkode_mri_get_current_time"

  external set_diagnostics : ('a, 'k) session -> Logfile.t -> unit
      = "sunml_arkode_mri_set_diagnostics"

  external clear_diagnostics : ('a, 'k) session -> unit
      = "sunml_arkode_mri_clear_diagnostics"

  external set_error_file : ('a, 'k) session -> Logfile.t -> unit
      = "sunml_arkode_mri_set_error_file"

  external c_set_err_handler_fn  : ('a, 'k) session -> unit
      = "sunml_arkode_mri_set_err_handler_fn"

  let set_err_handler_fn s ferrh =
    s.errh <- ferrh;
    c_set_err_handler_fn s

  external clear_err_handler_fn  : ('a, 'k) session -> unit
      = "sunml_arkode_mri_clear_err_handler_fn"

  let clear_err_handler_fn s =
    s.errh <- dummy_errh;
    clear_err_handler_fn s

  external set_fixed_step : ('d, 'k) session -> float -> unit
      = "sunml_arkode_mri_set_fixed_step"

  module Coupling = struct (* {{{ *)

    type cptr

    (* Synchronized with arkode_mri_coupling_index in arkode_ml.h *)
    type t = {
      cptr              : cptr;
      nmat              : int;
      stages            : int;
      method_order      : int;
      embedding_order   : int;
      coupling_matrices : RealArray.t array array;
      abscissae          : RealArray.t;
    }

    let nmat { nmat; _ } = nmat
    let stages { stages; _ } = stages
    let method_order { method_order; _ } = method_order
    let embedding_order { embedding_order; _ } = embedding_order
    let coupling_matrices { coupling_matrices; _ } = coupling_matrices
    let abscissae { abscissae; _ } = abscissae

    external c_make
      : int * int * int * int * RealArray.t array array * RealArray.t -> cptr
      = "sunml_arkode_mri_coupling_make"

    let make ~method_order ~embedding_order g c =
      let nmat = Array.length g in
      let stages = RealArray.length c in
      if nmat < 1 || stages < 1 then invalid_arg "zero-length array";
      let check ra = RealArray.length ra <> stages in
      if Array.exists (fun gi -> Array.length gi <> stages
                                 || Array.exists check gi) g
        then invalid_arg "coupling matrice incompatible with abscissae";
      {
        cptr = c_make (nmat, stages, method_order, embedding_order, g, c);
        nmat;
        stages;
        method_order;
        embedding_order;
        coupling_matrices = g;
        abscissae = c;
      }

    (* Synchronized with arkode_mri_coupling_table_tag in arkode_ml.h *)
    type coupling_table =
      | MIS_KW3
      | GARK_ERK45a
      | GARK_IRK21a
      | GARK_ESDIRK34a

    external load_table : coupling_table -> t
      = "sunml_arkode_mri_coupling_load_table"

    external c_mis_to_mri : int -> int -> ButcherTable.t -> t
      = "sunml_arkode_mri_coupling_mistomri"

    let mis_to_mri ~method_order ~embedding_order bt =
      c_mis_to_mri method_order embedding_order bt

    external copy : t -> t
      = "sunml_arkode_mri_coupling_copy"

    external space : t -> int * int
      = "sunml_arkode_mri_coupling_space"

    external write : t -> Logfile.t -> unit
      = "sunml_arkode_mri_coupling_write"

  end (* }}} *)

  external set_coupling : ('d, 'k) session -> Coupling.t -> unit
    = "sunml_arkode_mri_set_coupling"

  external c_set_table
    : ('d, 'k) session -> int -> ButcherTable.t option -> unit
    = "sunml_arkode_mri_set_table"

  let set_table s q bt =
    c_set_table s q (Some bt)

  external c_set_table_num : ('d, 'k) session -> int -> unit
      = "sunml_arkode_mri_set_table_num"

  let set_table_num s tn =
    c_set_table_num s (ButcherTable.int_of_erk_table tn)

  external set_defaults           : ('a, 'k) session -> unit
      = "sunml_arkode_mri_set_defaults"

  external set_interpolant_type : ('d, 'k) session -> interpolant_type -> unit
      = "sunml_arkode_mri_set_interpolant_type"

  external set_interpolant_degree : ('d, 'k) session -> int -> unit
      = "sunml_arkode_mri_set_interpolant_degree"

  external set_max_hnil_warns     : ('a, 'k) session -> int -> unit
      = "sunml_arkode_mri_set_max_hnil_warns"
  external set_max_num_steps      : ('a, 'k) session -> int -> unit
      = "sunml_arkode_mri_set_max_num_steps"
  external set_stop_time          : ('a, 'k) session -> float -> unit
      = "sunml_arkode_mri_set_stop_time"

  type 'd pre_inner_fn  = 'd Arkode_impl.Global.pre_inner_fn
  type 'd post_inner_fn = 'd Arkode_impl.Global.post_inner_fn

  external c_set_pre_inner_fn : ('a, 'k) session -> bool -> unit
      = "sunml_arkode_mri_set_pre_inner_fn"

  external c_set_post_inner_fn : ('a, 'k) session -> bool -> unit
      = "sunml_arkode_mri_set_post_inner_fn"

  let set_pre_inner_fn s f =
    s.preinnerfn <- f;
    c_set_pre_inner_fn s true

  let clear_pre_inner_fn s =
    c_set_pre_inner_fn s false;
    s.preinnerfn <- dummy_preinnerfn

  let set_post_inner_fn s f =
    s.postinnerfn <- f;
    c_set_post_inner_fn s true

  let clear_post_inner_fn s =
    c_set_post_inner_fn s false;
    s.postinnerfn <- dummy_postinnerfn

  external c_set_postprocess_step_fn : ('a, 'k) session -> bool -> unit
      = "sunml_arkode_mri_set_postprocess_step_fn"

  let set_postprocess_step_fn s fn =
    s.poststepfn <- fn;
    c_set_postprocess_step_fn s true

  let clear_postprocess_step_fn s =
    s.poststepfn <- dummy_poststepfn;
    c_set_postprocess_step_fn s false

  external set_nonlinear : ('a, 'k) session -> unit
    = "sunml_arkode_mri_set_nonlinear"

  external c_set_stage_predict_fn
      : ('d, 'k) session -> bool -> unit
      = "sunml_arkode_mri_set_stage_predict_fn"

  let set_stage_predict_fn s fn =
    s.stagepredictfn <- fn;
    c_set_stage_predict_fn s true

  let clear_stage_predict_fn s =
    c_set_stage_predict_fn s false;
    s.stagepredictfn <- dummy_stagepredictfn


  external c_set_root_direction   : ('a, 'k) session -> RootDirs.t -> unit
      = "sunml_arkode_mri_set_root_direction"

  let set_root_direction s rda =
    c_set_root_direction s (RootDirs.copy (get_num_roots s) rda)

  let set_all_root_directions s rd =
    c_set_root_direction s (RootDirs.make (get_num_roots s) rd)

  external set_no_inactive_root_warn      : ('a, 'k) session -> unit
      = "sunml_arkode_mri_set_no_inactive_root_warn"

  external get_current_state : ('d, 'k) session -> 'd
      = "sunml_arkode_mri_get_current_state"

  external get_current_coupling : ('d, 'k) session -> Coupling.t
    = "sunml_arkode_mri_get_current_coupling"

  external write_coupling : ('d, 'k) session -> Logfile.t -> unit
    = "sunml_arkode_mri_write_coupling"

  (* must correspond to arkode_nonlin_system_data_index in arkode_ml.h *)
  type 'd nonlin_system_data = {
    tcur  : float;
    zpred : 'd;
    zi    : 'd;
    fi    : 'd;
    gamma : float;
    sdata : 'd;
  }

  external get_nonlin_system_data
    : ('d, 'k) session -> 'd nonlin_system_data
    = "sunml_arkode_mri_get_nonlin_system_data"

  external compute_state
    : ('d, 'k) session
      -> ('d, 'k) Nvector.t
      -> ('d, 'k) Nvector.t
      -> unit
    = "sunml_arkode_mri_compute_state"

  external get_num_g_evals                : ('a, 'k) session -> int
      = "sunml_arkode_mri_get_num_g_evals"

  external write_parameters : ('d, 'k) session -> Logfile.t -> unit
      = "sunml_arkode_mri_write_parameters"

  external set_nonlin_conv_coef   : ('a, 'k) session -> float -> unit
      = "sunml_arkode_mri_set_nonlin_conv_coef"
  external set_nonlin_crdown      : ('a, 'k) session -> float -> unit
      = "sunml_arkode_mri_set_nonlin_crdown"
  external set_nonlin_rdiv        : ('a, 'k) session -> float -> unit
      = "sunml_arkode_mri_set_nonlin_rdiv"
  external set_delta_gamma_max    : ('a, 'k) session -> float -> unit
      = "sunml_arkode_mri_set_delta_gamma_max"
  external set_lsetup_frequency   : ('a, 'k) session -> int -> unit
      = "sunml_arkode_mri_set_lsetup_frequency"
  external set_predictor_method : ('d, 'k) session -> predictor_method -> unit
      = "sunml_arkode_mri_set_predictor_method"
  external set_max_nonlin_iters   : ('a, 'k) session -> int -> unit
      = "sunml_arkode_mri_set_max_nonlin_iters"

  external get_current_gamma : ('d, 'k) session -> float
      = "sunml_arkode_mri_get_current_gamma"

  external get_tol_scale_factor           : ('a, 'k) session -> float
      = "sunml_arkode_mri_get_tol_scale_factor"

  external c_get_err_weights
      : ('a, 'k) session -> ('a, 'k) nvector -> unit
      = "sunml_arkode_mri_get_err_weights"

  let get_err_weights s ew =
    if Sundials_configuration.safe then s.checkvec ew;
    c_get_err_weights s ew

  external get_num_lin_solv_setups : ('a, 'k) session -> int
      = "sunml_arkode_mri_get_num_lin_solv_setups"

  external get_num_nonlin_solv_iters      : ('a, 'k) session -> int
      = "sunml_arkode_mri_get_num_nonlin_solv_iters"

  external get_num_nonlin_solv_conv_fails : ('a, 'k) session -> int
      = "sunml_arkode_mri_get_num_nonlin_solv_conv_fails"

  external get_nonlin_solv_stats          : ('a, 'k) session -> int * int
      = "sunml_arkode_mri_get_nonlin_solv_stats"
end (* }}} *)

(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "sunml_arkode_init_module"

let _ =
  c_init_module
    (* Exceptions must be listed in the same order as
       arkode_exn_index.  *)
    [|IllInput;
      TooClose;
      TooMuchWork;
      TooMuchAccuracy;
      InnerStepFail None;
      ErrFailure;
      ConvergenceFailure;
      LinearInitFailure;
      LinearSetupFailure None;
      LinearSolveFailure None;
      NonlinearInitFailure;
      NonlinearSetupFailure;
      NonlinearSetupRecoverable;
      NonlinearOperationError;
      MassInitFailure;
      MassSetupFailure None;
      MassSolveFailure None;
      MassMultFailure;
      RhsFuncFailure;
      FirstRhsFuncFailure;
      RepeatedRhsFuncFailure;
      UnrecoverableRhsFuncFailure;
      RootFuncFailure;
      PostprocStepFailure;
      BadK;
      BadT;
      VectorOpErr;
    |]
