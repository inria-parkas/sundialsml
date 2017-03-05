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

include Kinsol_impl
open Sundials

exception IllInput                       (* KIN_ILL_INPUT *)
exception LineSearchNonConvergence       (* KIN_LINESEARCH_NONCONV *)
exception MaxIterationsReached           (* KIN_MAXITER_REACHED *)
exception MaxNewtonStepExceeded          (* KIN_MXNEWT_5X_EXCEEDED *)
exception LineSearchBetaConditionFailure (* KIN_LINESEARCH_BCFAIL *)
exception LinearSolverNoRecovery         (* KIN_LINSOLV_NO_RECOVERY *)
exception LinearSolverInitFailure        (* KIN_LINIT_FAIL *)
exception LinearSetupFailure             (* KIN_LSETUP_FAIL *)
exception LinearSolverFailure            (* KIN_LSOLVE_FAIL *)
exception SystemFunctionFailure          (* KIN_SYSFUNC_FAIL *)
exception FirstSystemFunctionFailure     (* KIN_FIRST_SYSFUNC_FAIL *)
exception RepeatedSystemFunctionFailure  (* KIN_REPTD_SYSFUNC_ERR *)
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

type 'k serial_linear_solver = (RealArray.t, 'k) linear_solver
                               constraint 'k = [>Nvector_serial.kind]
type 'a double = 'a * 'a

(* interface *)

let int_default = function None -> 0 | Some v -> v

module Dls =
  struct
    include DlsTypes

    external c_dls_dense : 'k serial_session -> bool -> unit
      = "c_kinsol_dls_dense"

    external c_dls_lapack_dense : 'k serial_session -> bool -> unit
      = "c_kinsol_dls_lapack_dense"

    external c_dls_band : 'k serial_session -> int -> int -> bool -> unit
      = "c_kinsol_dls_band"

    external c_dls_lapack_band : 'k serial_session -> int -> int -> bool -> unit
      = "c_kinsol_dls_lapack_band"

    let dense ?jac () s nv =
      s.neqs <- Sundials.RealArray.length (Nvector.unwrap nv);
      c_dls_dense s (jac <> None);
      s.ls_precfns <- NoPrecFns;
      s.ls_callbacks <- match jac with
                        | None   -> DlsDenseCallback no_dense_callback
                        | Some f -> DlsDenseCallback { jacfn = f; dmat = None }

    let lapack_dense ?jac () s nv =
      s.neqs <- Sundials.RealArray.length (Nvector.unwrap nv);
      c_dls_lapack_dense s (jac <> None);
      s.ls_precfns <- NoPrecFns;
      s.ls_callbacks <- match jac with
                        | None   -> DlsDenseCallback no_dense_callback
                        | Some f -> DlsDenseCallback { jacfn = f; dmat = None }

    let band ?jac { mupper; mlower } s nv =
      s.neqs <- Sundials.RealArray.length (Nvector.unwrap nv);
      c_dls_band s mupper mlower (jac <> None);
      s.ls_precfns <- NoPrecFns;
      s.ls_callbacks <- match jac with
                        | None   -> DlsBandCallback no_band_callback
                        | Some f -> DlsBandCallback { bjacfn = f; bmat = None }

    let lapack_band ?jac { mupper; mlower } s nv =
      s.neqs <- Sundials.RealArray.length (Nvector.unwrap nv);
      c_dls_lapack_band s mupper mlower (jac <> None);
      s.ls_precfns <- NoPrecFns;
      s.ls_callbacks <- match jac with
                        | None   -> DlsBandCallback no_band_callback
                        | Some f -> DlsBandCallback { bjacfn = f; bmat = None }

    let invalidate_callback session =
      match session.ls_callbacks with
      | DlsDenseCallback ({ dmat = Some d } as cb) ->
          Dls.DenseMatrix.invalidate d;
          cb.dmat <- None
      | DlsBandCallback  ({ bmat = Some d } as cb) ->
          Dls.BandMatrix.invalidate d;
          cb.bmat <- None
      | _ -> ()

    external set_dense_jac_fn : 'k serial_session -> unit
        = "c_kinsol_dls_set_dense_jac_fn"

    let set_dense_jac_fn s fjacfn =
      match s.ls_callbacks with
      | DlsDenseCallback _ ->
          invalidate_callback s;
          s.ls_callbacks <- DlsDenseCallback { jacfn = fjacfn; dmat = None };
          set_dense_jac_fn s
      | _ -> raise Sundials.InvalidLinearSolver

    external clear_dense_jac_fn : 'k serial_session -> unit
        = "c_kinsol_dls_clear_dense_jac_fn"

    let clear_dense_jac_fn s =
      match s.ls_callbacks with
      | DlsDenseCallback _ ->
          invalidate_callback s;
          s.ls_callbacks <- DlsDenseCallback no_dense_callback;
          clear_dense_jac_fn s
      | _ -> raise Sundials.InvalidLinearSolver

    external set_band_jac_fn : 'k serial_session -> unit
        = "c_kinsol_dls_set_band_jac_fn"

    let set_band_jac_fn s fbandjacfn =
      match s.ls_callbacks with
      | DlsBandCallback _ ->
          invalidate_callback s;
          s.ls_callbacks <-
            DlsBandCallback { bjacfn = fbandjacfn; bmat = None };
          set_band_jac_fn s
      | _ -> raise Sundials.InvalidLinearSolver

    external clear_band_jac_fn : 'k serial_session -> unit
        = "c_kinsol_dls_clear_band_jac_fn"

    let clear_band_jac_fn s =
      match s.ls_callbacks with
      | DlsBandCallback _ ->
          invalidate_callback s;
          s.ls_callbacks <- DlsBandCallback no_band_callback;
          clear_band_jac_fn s
      | _ -> raise Sundials.InvalidLinearSolver

    external get_work_space : 'k serial_session -> int * int
        = "c_kinsol_dls_get_work_space"

    let get_work_space s =
      ls_check_dls s;
      get_work_space s

    external get_num_jac_evals : 'k serial_session -> int
        = "c_kinsol_dls_get_num_jac_evals"

    let get_num_jac_evals s =
      ls_check_dls s;
      get_num_jac_evals s

    external get_num_func_evals : 'k serial_session -> int
        = "c_kinsol_dls_get_num_func_evals"

    let get_num_func_evals s =
      ls_check_dls s;
      get_num_func_evals s
  end

module Sls =
  struct
    include SlsTypes

    module Klu = struct

      (* Must correspond with kin_klu_ordering_tag *)
      type ordering =
           Amd
         | ColAmd
         | Natural

      external c_klu
        : 'k serial_session -> Sls_impl.sformat -> int -> int -> unit
        = "c_kinsol_klu_init"

      let solver sformat f nnz session nv =
        if not Sundials_config.klu_enabled
          then raise Sundials.NotImplementedBySundialsVersion;
        session.neqs <- Sundials.RealArray.length (Nvector.unwrap nv);
        session.ls_callbacks <- SlsKluCallback { jacfn = f; smat = None };
        session.ls_precfns <- NoPrecFns;
        c_klu session sformat session.neqs nnz

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
        = "c_kinsol_klu_set_ordering"

      let set_ordering session ordering =
        ls_check_klu session;
        c_set_ordering session ordering

      external c_reinit : 'k serial_session -> int -> int -> bool -> unit
        = "c_kinsol_klu_reinit"

      let reinit session n nnz realloc =
        ls_check_klu session;
        c_reinit session n nnz realloc

      external c_get_num_jac_evals : 'k serial_session -> int
        = "c_kinsol_klu_get_num_jac_evals"

      let get_num_jac_evals session =
        ls_check_klu session;
        c_get_num_jac_evals session

    end

    module Superlumt = struct

      (* Must correspond with kinsol_superlumt_ordering_tag *)
      type ordering =
           Natural
         | MinDegreeProd
         | MinDegreeSum
         | ColAmd

      external c_superlumt : 'k serial_session -> int -> int -> int -> unit
        = "c_kinsol_superlumt_init"

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
        = "c_kinsol_superlumt_set_ordering"

      let set_ordering session ordering =
        ls_check_superlumt session;
        c_set_ordering session ordering

      external c_get_num_jac_evals : 'k serial_session -> int
        = "c_kinsol_superlumt_get_num_jac_evals"

      let get_num_jac_evals session =
        ls_check_superlumt session;
        c_get_num_jac_evals session

    end
  end

module Spils =
  struct
    include SpilsTypes

    external c_spgmr : ('a, 'k) session -> int -> unit
      = "c_kinsol_spils_spgmr"

    external c_spfgmr : ('a, 'k) session -> int -> unit
      = "c_kinsol_spils_spfgmr"

    external c_spbcg : ('a, 'k) session -> int -> unit
      = "c_kinsol_spils_spbcg"

    external c_sptfqmr : ('a, 'k) session -> int -> unit
      = "c_kinsol_spils_sptfqmr"

    external c_set_max_restarts     : ('a, 'k) session -> int -> unit
        = "c_kinsol_spils_set_max_restarts"

    external c_set_preconditioner : ('a, 'k) session -> bool -> bool -> unit
        = "c_kinsol_spils_set_preconditioner"

    external c_set_jac_times_vec_fn : ('a, 'k) session -> bool -> unit
        = "c_kinsol_spils_set_jac_times_vec_fn"

    let init_preconditioner solve setup session nv =
      c_set_preconditioner session (solve <> None) (setup <> None);
      session.ls_precfns <- PrecFns { prec_solve_fn = solve;
                                      prec_setup_fn = setup }

    let prec_none = InternalPrecNone (fun session nv ->
        session.ls_precfns <- NoPrecFns)
    let prec_right ?setup ?solve () =
      InternalPrecRight (init_preconditioner solve setup)

    let set_jac_times_vec_fn s f =
      match s.ls_callbacks with
      | SpilsCallback _ ->
          s.ls_callbacks <- SpilsCallback (Some f);
          c_set_jac_times_vec_fn s true
      | _ -> raise Sundials.InvalidLinearSolver

    let init_spils init maxl jac_times_vec prec session nv =
      init session maxl;
      (match prec with
       | InternalPrecNone set_prec  -> set_prec session nv
       | InternalPrecRight set_prec -> set_prec session nv);
      session.ls_callbacks <- SpilsCallback jac_times_vec;
      (match jac_times_vec with
       | None -> ()
       | Some jtv -> set_jac_times_vec_fn session jtv)

    let spgmr ?(maxl=0) ?(maxr=5) ?jtv prec session nv =
      init_spils c_spgmr maxl jtv prec session nv;
      (* Note: we can skip set_max_restarts only when initializing a
         fresh solver.  *)
      if maxr <> 5 then c_set_max_restarts session maxr

    let spfgmr ?(maxl=0) ?(maxr=5) ?jtv prec session nv =
      init_spils c_spfgmr maxl jtv prec session nv;
      (* Note: we can skip set_max_restarts only when initializing a
         fresh solver.  *)
      if maxr <> 5 then c_set_max_restarts session maxr

    let spbcg ?(maxl=0) ?jtv prec session nv =
      init_spils c_spbcg maxl jtv prec session nv

    let sptfqmr ?(maxl=0) ?jtv prec session nv =
      init_spils c_sptfqmr maxl jtv prec session nv

    let set_preconditioner s ?setup ?solve () =
      match s.ls_callbacks with
      | SpilsCallback _ ->
          c_set_preconditioner s (solve <> None) (setup <> None);
          s.ls_precfns <- PrecFns { prec_setup_fn = setup;
                                    prec_solve_fn = solve }
      | _ -> raise Sundials.InvalidLinearSolver

    let clear_jac_times_vec_fn s =
      match s.ls_callbacks with
      | SpilsCallback cbs ->
          c_set_jac_times_vec_fn s false;
          s.ls_callbacks <- SpilsCallback None
      | _ -> raise Sundials.InvalidLinearSolver

    external get_work_space       : ('a, 'k) session -> int * int
        = "c_kinsol_spils_get_work_space"

    let get_work_space s =
      ls_check_spils s;
      get_work_space s

    external get_num_lin_iters    : ('a, 'k) session -> int
        = "c_kinsol_spils_get_num_lin_iters"

    let get_num_lin_iters s =
      ls_check_spils s;
      get_num_lin_iters s

    external get_num_conv_fails   : ('a, 'k) session -> int
        = "c_kinsol_spils_get_num_conv_fails"

    let get_num_conv_fails s =
      ls_check_spils s;
      get_num_conv_fails s

    external get_num_prec_evals   : ('a, 'k) session -> int
        = "c_kinsol_spils_get_num_prec_evals"

    let get_num_prec_evals s =
      ls_check_spils s;
      get_num_prec_evals s

    external get_num_prec_solves  : ('a, 'k) session -> int
        = "c_kinsol_spils_get_num_prec_solves"

    let get_num_prec_solves s =
      ls_check_spils s;
      get_num_prec_solves s

    external get_num_jtimes_evals : ('a, 'k) session -> int
        = "c_kinsol_spils_get_num_jtimes_evals"

    let get_num_jtimes_evals s =
      ls_check_spils s;
      get_num_jtimes_evals s

    external get_num_func_evals    : ('a, 'k) session -> int
        = "c_kinsol_spils_get_num_func_evals"

    let get_num_func_evals s =
      ls_check_spils s;
      get_num_func_evals s
  end

module Alternate =
  struct
    include AlternateTypes

    external c_set_alternate
      : ('data, 'kind) session -> bool -> bool -> unit
      = "c_kinsol_set_alternate"

    external get_u_uscale : ('data, 'kind) session -> 'data * 'data
      = "c_kinsol_get_u_uscale"

    external get_f_fscale  : ('data, 'kind) session -> 'data * 'data
      = "c_kinsol_get_f_fscale"

    external set_sjpnorm   : ('data, 'kind) session -> float -> unit
      = "c_kinsol_set_sjpnorm"

    external set_sfdotjp   : ('data, 'kind) session -> float -> unit
      = "c_kinsol_set_sfdotjp"

    let make f s nv =
      let { linit; lsetup; lsolve } as cb = f s nv in
      c_set_alternate s (linit <> None) (lsetup <> None);
      s.ls_precfns <- NoPrecFns;
      s.ls_callbacks <- AlternateCallback cb

  end

external set_error_file : ('a, 'k) session -> Sundials.Logfile.t -> unit
    = "c_kinsol_set_error_file"

external c_set_err_handler_fn : ('a, 'k) session -> unit
    = "c_kinsol_set_err_handler_fn"

let set_err_handler_fn s ferrh =
  s.errh <- ferrh;
  c_set_err_handler_fn s

external c_clear_err_handler_fn : ('a, 'k) session -> unit
    = "c_kinsol_clear_err_handler_fn"

let clear_err_handler_fn s =
  s.errh <- dummy_errh;
  c_clear_err_handler_fn s

external set_info_file : ('a, 'k) session -> Sundials.Logfile.t -> unit
    = "c_kinsol_set_info_file"

external c_set_info_handler_fn : ('a, 'k) session -> unit
    = "c_kinsol_set_info_handler_fn"

let set_info_handler_fn s finfoh =
  s.infoh <- finfoh;
  c_set_info_handler_fn s

external c_clear_info_handler_fn : ('a, 'k) session -> unit
    = "c_kinsol_clear_info_handler_fn"

let clear_info_handler_fn s =
  s.infoh <- dummy_infoh;
  c_clear_info_handler_fn s

external set_print_level : ('a, 'k) session -> print_level -> unit
    = "c_kinsol_set_print_level"

external set_num_max_iters : ('a, 'k) session -> int -> unit
    = "c_kinsol_set_num_max_iters"

external set_maa : ('a, 'k) session -> int -> unit
    = "c_kinsol_set_maa"

external c_set_no_init_setup : ('a, 'k) session -> bool -> unit
    = "c_kinsol_set_no_init_setup"

let set_no_init_setup s = c_set_no_init_setup s true
let set_init_setup s = c_set_no_init_setup s false

external c_set_no_res_mon : ('a, 'k) session -> bool -> unit
    = "c_kinsol_set_no_res_mon"

let set_no_res_mon s = c_set_no_res_mon s true
let set_res_mon s = c_set_no_res_mon s false

external set_max_sub_setup_calls : ('a, 'k) session -> int -> unit
    = "c_kinsol_set_max_sub_setup_calls"

external set_max_setup_calls : ('a, 'k) session -> int -> unit
    = "c_kinsol_set_max_setup_calls"

external c_set_eta_form : ('a, 'k) session -> eta_choice -> unit
    = "c_kinsol_set_eta_form"

external c_set_eta_const_value : ('a, 'k) session -> float -> unit
    = "c_kinsol_set_eta_form"

external c_set_eta_params : ('a, 'k) session -> float -> float -> unit
    = "c_kinsol_set_eta_form"

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
    = "c_kinsol_set_res_mon_const_value"

external c_set_res_mon_params : ('a, 'k) session -> float -> float -> unit
    = "c_kinsol_set_res_mon_params"

let set_res_mon_params s ?omegamin ?omegamax () =
  c_set_res_mon_params s (float_default omegamin) (float_default omegamax)

external c_set_no_min_eps : ('a, 'k) session -> bool -> unit
    = "c_kinsol_set_no_min_eps"

let set_no_min_eps s = c_set_no_min_eps s true
let set_min_eps s = c_set_no_min_eps s false

external set_max_newton_step : ('a, 'k) session -> float -> unit
    = "c_kinsol_set_max_newton_step"

external set_max_beta_fails : ('a, 'k) session -> float -> unit
    = "c_kinsol_set_max_beta_fails"

external set_rel_err_func : ('a, 'k) session -> float -> unit
    = "c_kinsol_set_rel_err_func"

external set_func_norm_tol : ('a, 'k) session -> float -> unit
    = "c_kinsol_set_func_norm_tol"

external set_scaled_step_tol : ('a, 'k) session -> float -> unit
    = "c_kinsol_set_scaled_step_tol"

external c_set_constraints : ('a, 'k) session -> ('a, 'k) nvector -> unit
    = "c_kinsol_set_constraints"

let set_constraints s cc =
  if Sundials_config.safe then s.checkvec cc;
  c_set_constraints s cc

let set_sys_func s fsys =
  s.sysfn <- fsys

external get_work_space : ('a, 'k) session -> int * int
    = "c_kinsol_get_work_space"

external get_num_func_evals : ('a, 'k) session -> int
    = "c_kinsol_get_num_func_evals"

external get_num_nonlin_solv_iters : ('a, 'k) session -> int
    = "c_kinsol_get_num_nonlin_solv_iters"

external get_num_beta_cond_fails : ('a, 'k) session -> int
    = "c_kinsol_get_num_beta_cond_fails"

external get_num_backtrack_ops : ('a, 'k) session -> int
    = "c_kinsol_get_num_backtrack_ops"

external get_func_norm : ('a, 'k) session -> float
    = "c_kinsol_get_func_norm"

external get_step_length : ('a, 'k) session -> float
    = "c_kinsol_get_step_length"

external c_init
    : ('a, 'k) session Weak.t -> ('a, 'k) nvector -> int option -> int option
      -> (kin_mem * c_weak_ref)
    = "c_kinsol_init"

external c_session_finalize : ('a, 'k) session -> unit
    = "c_kinsol_session_finalize"

let session_finalize s =
  Dls.invalidate_callback s;
  c_session_finalize s

let init ?maxi ?maa ?lsolver f u0 =
  let checkvec = Nvector.check u0 in
  let weakref = Weak.create 1 in
  let kin_mem, backref = c_init weakref u0 maxi maa
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

          ls_callbacks = NoCallbacks;
          ls_precfns   = NoPrecFns;
        } in
  Gc.finalise session_finalize session;
  Weak.set weakref 0 (Some session);
  (match lsolver with Some lsolver -> lsolver session u0 | None -> ());
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
    = "c_kinsol_solve"

let solve s u strategy u_scale f_scale =
  if Sundials_config.safe then
    (s.checkvec u;
     s.checkvec u_scale;
     s.checkvec f_scale);
  if strategy <> FixedPoint && s.ls_callbacks = NoCallbacks
  then raise MissingLinearSolver
  else c_solve s u strategy u_scale f_scale

(* Let C code know about some of the values in this module.  *)
external c_init_module : 'fcns -> exn array -> unit =
  "c_kinsol_init_module"

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
      LinearSetupFailure;
      LinearSolverFailure;
      SystemFunctionFailure;
      FirstSystemFunctionFailure;
      RepeatedSystemFunctionFailure;
    |]
