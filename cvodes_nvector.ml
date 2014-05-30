(***********************************************************************)
(*                                                                     *)
(*               OCaml interface to (serial) Sundials                  *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

include Cvode_session_nvector

let add_fwdsensext s =
  match s.sensext with
  | FwdSensExt se -> ()
  | BwdSensExt _ -> failwith "Quadrature.add_fwdsensext: internal error"
  | NoSensExt ->
      s.sensext <- FwdSensExt {
        quadrhsfn     = (fun _ _ _ -> ());
        senspvals     = None;
        sensrhsfn     = (fun _ _ _ _ _ _ _ -> ());
        sensrhsfn1    = (fun _ _ _ _ _ _ _ _ -> ());
        quadsensrhsfn = (fun _ _ _ _ _ _ _ -> ());
      }

(* TODO: add callback 'trampolines' *)

let _ = List.iter (fun (nm, ex) -> Callback.register_exception nm ex)
  [
    ("cvodes_RecoverableFailure",      Sundials.RecoverableFailure);
  ]

module Quadrature =
  struct
    exception QuadNotInitialized
    exception QuadRhsFuncFailure
    exception FirstQuadRhsFuncErr
    exception RepeatedQuadRhsFuncErr
    exception UnrecoverableQuadRhsFuncErr

    let _ = List.iter (fun (nm, ex) -> Callback.register_exception nm ex)
      [
        ("cvodes_QuadNotInitialized",           QuadNotInitialized);
        ("cvodes_QuadRhsFuncFailure",           QuadRhsFuncFailure);
        ("cvodes_FirstQuadRhsFuncErr",          FirstQuadRhsFuncErr);
        ("cvodes_RepeatedQuadRhsFuncErr",       RepeatedQuadRhsFuncErr);
        ("cvodes_UnrecoverableQuadRhsFuncErr",  UnrecoverableQuadRhsFuncErr)
      ]

    let fwdsensext s =
      match s.sensext with
      | FwdSensExt se -> se
      | _ -> raise QuadNotInitialized

    type 'a quadrhsfn = float -> 'a -> 'a -> unit

    external c_init : 'a session -> 'a nvector -> unit
        = "c_nvec_cvodes_quad_init"

    let init s f v0 =
      add_fwdsensext s;
      let se = fwdsensext s in
      se.quadrhsfn <- f;
      c_init s v0

    external reinit : 'a session -> 'a nvector -> unit
      = "c_nvec_cvodes_quad_reinit"

    external set_err_con    : 'a session -> bool -> unit
        = "c_nvec_cvodes_quad_set_err_con"
    external sv_tolerances  : 'a session -> float -> 'a nvector -> unit
        = "c_nvec_cvodes_quad_sv_tolerances"
    external ss_tolerances  : 'a session -> float -> float -> unit
        = "c_cvodes_quad_ss_tolerances"

    type 'a tolerance =
        NoStepSizeControl
      | SSTolerances of float * float
      | SVTolerances of float * 'a nvector

    let set_tolerances s tol =
      match tol with
      | NoStepSizeControl -> set_err_con s false
      | SSTolerances (rel, abs) -> (set_err_con s true;
                                    ss_tolerances s rel abs)
      | SVTolerances (rel, abs) -> (set_err_con s true;
                                    sv_tolerances s rel abs)

    external get : 'a session -> 'a nvector -> float
        = "c_cvodes_quad_get"

    external get_dky : 'a session -> float -> int -> 'a nvector -> unit
        = "c_cvodes_quad_get_dky"

    external get_num_rhs_evals       : 'a session -> int
        = "c_cvodes_quad_get_num_rhs_evals"

    external get_num_err_test_fails  : 'a session -> int
        = "c_cvodes_quad_get_num_err_test_fails"

    external get_quad_err_weights : 'a session -> 'a nvector -> unit
        = "c_cvodes_quad_get_quad_err_weights"

    external get_stats : 'a session -> int * int
        = "c_cvodes_quad_get_stats"
  end

module Sensitivity =
  struct
    type 'a tolerance =
        SSTolerances of float * float
      | SVTolerances of float * 'a nvector
      | EETolerances

    external set_err_con : 'a session -> bool -> unit
        = "c_cvodes_sens_set_err_con"

    external ss_tolerances  : 'a session -> float -> float -> unit
        = "c_cvodes_sens_ss_tolerances"

    external ee_tolerances  : 'a session -> unit
        = "c_nvec_cvodes_sens_ee_tolerances"

    external sv_tolerances  : 'a session -> float -> 'a nvector -> unit
        = "c_nvec_cvodes_sens_sv_tolerances"

    let set_tolerances s tol =
      match tol with
      | SSTolerances (rel, abs) -> ss_tolerances s rel abs
      | SVTolerances (rel, abs) -> sv_tolerances s rel abs
      | EETolerances -> ee_tolerances s
    

    exception SensNotInitialized
    exception SensRhsFuncFailure
    exception FirstSensRhsFuncErr
    exception RepeatedSensRhsFuncErr
    exception UnrecoverableSensRhsFuncErr
    exception BadIS

    let _ = List.iter (fun (nm, ex) -> Callback.register_exception nm ex)
      [
        ("cvodes_SensNotInitialized",           SensNotInitialized);
        ("cvodes_SensRhsFuncFailure",           SensRhsFuncFailure);
        ("cvodes_FirstSensRhsFuncErr",          FirstSensRhsFuncErr);
        ("cvodes_RepeatedSensRhsFuncErr",       RepeatedSensRhsFuncErr);
        ("cvodes_UnrecoverableSensRhsFuncErr",  UnrecoverableSensRhsFuncErr);
        ("cvodes_BadIS",                        BadIS)
      ]

    let fwdsensext s =
      match s.sensext with
      | FwdSensExt se -> se
      | _ -> raise SensNotInitialized

    type sens_method =
        Simultaneous
      | Staggered
      | Staggered1

    type 'a _sensrhsfn = 'a sensrhsfn =
        AllAtOnce of
          (float -> 'a -> 'a -> 'a array -> 'a array -> 'a -> 'a -> unit)
      | OneByOne of
          (float -> 'a -> 'a -> int -> 'a -> 'a -> 'a -> 'a -> unit)
    type 'a sensrhsfn = 'a _sensrhsfn =
        AllAtOnce of
          (float -> 'a -> 'a -> 'a array -> 'a array -> 'a -> 'a -> unit)
      | OneByOne of
          (float -> 'a -> 'a -> int -> 'a -> 'a -> 'a -> 'a -> unit)

    type sens_params = {
        pvals  : Sundials.real_array option;
        pbar   : Sundials.real_array option;
        plist  : int array option;
      }

    let no_sens_params = { pvals = None; pbar = None; plist = None }

    external c_init : 'a session -> sens_method -> 'a nvector array -> unit
        = "c_nvec_cvodes_sens_init"

    external c_init_1 : 'a session -> sens_method -> 'a nvector array -> unit
        = "c_nvec_cvodes_sens_init_1"

    (* TODO: check that pbar and plist are ns long. *)
    external set_params : 'a session -> sens_params -> unit
        = "c_cvodes_sens_set_params"

    let init s tol fmethod sparams fm v0 =
      add_fwdsensext s;
      let se = fwdsensext s in
      c_init s fmethod v0;
      (match fm with
       | AllAtOnce f -> begin
           if fmethod = Staggered1 then
             failwith "Forward.init: Cannot combine AllAtOnce and Staggered1";
           se.sensrhsfn <- f;
           c_init s fmethod v0
         end
       | OneByOne f -> begin
           se.sensrhsfn1 <- f;
           c_init_1 s fmethod v0
        end);
      se.senspvals <- sparams.pvals;
      set_params s sparams;
      set_tolerances s tol

    external reinit : 'a session -> sens_method -> 'a nvector array -> unit
        = "c_nvec_cvodes_sens_reinit"

    external toggle_off : 'a session -> unit
        = "c_cvodes_sens_toggle_off"

    external get : 'a session -> 'a nvector array -> float
        = "c_nvec_cvodes_sens_get"

    external get_dky : 'a session -> float -> int -> 'a nvector array -> unit
        = "c_nvec_cvodes_sens_get_dky"

    external get' : 'a session -> int -> 'a nvector -> float
        = "c_nvec_cvodes_sens_get1"

    external get_dky' : 'a session -> float -> int -> int -> 'a nvector -> unit
        = "c_nvec_cvodes_sens_get_dky1"

    type dq_method = DQCentered | DQForward

    external set_dq_method : 'a session -> dq_method -> float -> unit
        = "c_cvodes_sens_set_dq_method"

    external set_max_nonlin_iters : 'a session -> int -> unit
        = "c_cvodes_sens_set_max_nonlin_iters"

    external get_num_sens_evals : 'a session -> int
        = "c_cvodes_sens_get_num_sens_evals"

    external get_num_rhs_evals : 'a session -> int
        = "c_cvodes_sens_get_num_rhs_evals"

    external get_num_err_test_fails : 'a session -> int
        = "c_cvodes_sens_get_num_err_test_fails"

    external get_num_lin_solv_setups : 'a session -> int
        = "c_cvodes_sens_get_num_lin_solv_setups"

    type sensitivity_stats = {
        num_rhs_evals : int;
        num_sens_evals :int;
        num_err_test_fails : int;
        num_lin_solv_setups :int;
      }

    external get_stats : 'a session -> sensitivity_stats
        = "c_cvodes_sens_get_stats"

    external get_err_weights : 'a session -> 'a nvector array -> unit
        = "c_nvec_cvodes_sens_get_err_weights"

    external get_num_nonlin_solv_iters : 'a session -> int
        = "c_cvodes_sens_get_num_nonlin_solv_iters"

    external get_num_nonlin_solv_conv_fails : 'a session -> int
        = "c_cvodes_sens_get_num_nonlin_solv_conf_fails"

    type nonlin_stats = {
        num_nonlin_solv_iters : int;
        num_nonlin_solv_conv_fails : int;
      }

    external get_nonlin_solv_stats : 'a session -> nonlin_stats
        = "c_cvodes_sens_get_nonlin_solv_stats"

    external get_num_stgr_nonlin_solv_iters
        : 'a session -> Sundials.lint_array -> unit
        = "c_cvodes_sens_get_num_stgr_nonlin_solv_iters"

    external get_num_stgr_nonlin_solv_conv_fails
        : 'a session -> Sundials.lint_array -> unit
        = "c_cvodes_sens_get_num_stgr_nonlin_solv_conv_fails"

    module Quadrature =
      struct

        exception QuadSensNotInitialized
        exception QuadSensRhsFuncFailure
        exception FirstQuadSensRhsFuncErr
        exception RepeatedQuadSensRhsFuncErr
        exception UnrecoverableQuadSensRhsFuncErr

        let _ = List.iter (fun (nm, ex) -> Callback.register_exception nm ex)
          [
            ("cvodes_QuadSensNotInitialized",     QuadSensNotInitialized);
            ("cvodes_QuadSensRhsFuncFailure",     QuadSensRhsFuncFailure);
            ("cvodes_FirstQuadSensRhsFuncErr",    FirstQuadSensRhsFuncErr);
            ("cvodes_RepeatedQuadSensRhsFuncErr", RepeatedQuadSensRhsFuncErr);
            ("cvodes_UnrecoverableQuadSensRhsFuncErr",
                                             UnrecoverableQuadSensRhsFuncErr)
          ]

        type 'a quadsensrhsfn =
           float                  (* t *)
           -> 'a nvector          (* y *)
           -> 'a nvector          (* yS *)
           -> 'a nvector          (* yQdot *)
           -> 'a nvector array    (* rhsvalQs *)
           -> 'a nvector          (* tmp1 *)
           -> 'a nvector          (* tmp2 *)
           -> unit

        external c_init : 'a session -> 'a nvector array -> unit
            = "c_nvec_cvodes_quadsens_init"

        let init s f v0 =
          let se = fwdsensext s in
          se.quadsensrhsfn <- f;
          c_init s v0

        external reinit : 'a session -> 'a nvector array -> unit
            = "c_nvec_cvodes_quadsens_reinit"

        type 'a tolerance =
            NoStepSizeControl
          | SSTolerances of float * float
          | SVTolerances of float * 'a nvector
          | EETolerances

        external set_err_con : 'a session -> bool -> unit
            = "c_cvodes_quadsens_set_err_con"

        external ss_tolerances  : 'a session -> float -> float -> unit
            = "c_cvodes_quadsens_ss_tolerances"

        external sv_tolerances  : 'a session -> float -> 'a nvector -> unit
            = "c_nvec_cvodes_quadsens_sv_tolerances"

        external ee_tolerances  : 'a session -> unit
            = "c_nvec_cvodes_quadsens_ee_tolerances"

        let set_tolerances s tol =
          match tol with
          | NoStepSizeControl -> set_err_con s false
          | SSTolerances (rel, abs) -> (set_err_con s true;
                                        ss_tolerances s rel abs)
          | SVTolerances (rel, abs) -> (set_err_con s true;
                                        sv_tolerances s rel abs)
          | EETolerances -> (set_err_con s true;
                             ee_tolerances s)
    
        external get : 'a session -> 'a nvector array -> float
            = "c_nvec_cvodes_quadsens_get"

        external get1 : 'a session -> int -> 'a nvector -> float
            = "c_nvec_cvodes_quadsens_get1"

        external get_dky : 'a session -> float -> int -> 'a nvector array -> unit
            = "c_nvec_cvodes_quadsens_get_dky"

        external get_dky1 : 'a session -> float -> int -> int -> 'a nvector array -> unit
            = "c_nvec_cvodes_quadsens_get_dky1"

        external get_num_rhs_evals       : 'a session -> int
            = "c_cvodes_quadsens_get_num_rhs_evals"

        external get_num_err_test_fails  : 'a session -> int
            = "c_cvodes_quadsens_get_num_err_test_fails"

        external get_quad_err_weights : 'a session -> 'a nvector array -> unit
            = "c_cvodes_quadsens_get_err_weights"

        external get_stats : 'a session -> int * int
            = "c_cvodes_quadsens_get_stats"
      end
    end

module Adjoint =
  struct
    exception AdjointNotInitialized
    exception NoForwardCall
    exception ForwardReinitializationFailed
    exception ForwardFailed
    exception NoBackwardProblem
    exception BadTB0
    exception BadT

    let _ = List.iter (fun (nm, ex) -> Callback.register_exception nm ex)
      [
        ("cvodes_AdjointNotInitialized",         AdjointNotInitialized);
        ("cvodes_NoForwardCall",                 NoForwardCall);
        ("cvodes_ForwardReinitializationFailed", ForwardReinitializationFailed);
        ("cvodes_ForwardFailed",                 ForwardFailed);
        ("cvodes_NoBackwardProblem",             NoBackwardProblem);
        ("cvodes_BadTB0",                        BadTB0);
        ("cvodes_BadT",                          BadT)
      ]

    type interpolation = IPolynomial | IHermite

    external init : 'a session -> int -> interpolation -> unit
        = "c_cvodes_adj_init"

    external forward_normal : 'a session -> float -> 'a nvector -> float * int
        = "c_nvec_cvodes_adj_forward_normal"

    external forward_one_step : 'a session -> float -> 'a nvector -> float * int
        = "c_nvec_cvodes_adj_forward_one_step"

    type 'a _brhsfn = 'a brhsfn =
        BackBasic of (float -> 'a -> 'a -> 'a -> unit)
      | BackWithSens of (float -> 'a -> 'a array -> 'a -> 'a -> unit)
    type 'a brhsfn = 'a _brhsfn =
        BackBasic of (float -> 'a -> 'a -> 'a -> unit)
      | BackWithSens of (float -> 'a -> 'a array -> 'a -> 'a -> unit)

    type 'a single_tmp = 'a
    type 'a triple_tmp = 'a * 'a * 'a

    type ('t, 'a) jacobian_arg = ('t, 'a) bjacobian_arg =
      {
        jac_t   : float;
        jac_y   : 'a;
        jac_yb  : 'a;
        jac_fyb : 'a;
        jac_tmp : 't
      }

    type 'a iter =
      | Newton of 'a linear_solver
      | Functional

    and 'a linear_solver =
(*
      (* TODO: not for nvectors *)
      | Dense of dense_jac_fnb option
      (* TODO: not for nvectors *)
      | LapackDense of dense_jac_fnb option
      (* TODO: not for nvectors *)
      | Band of bandrange * band_jac_fnb option
      (* TODO: not for nvectors *)
      | LapackBand of bandrange * band_jac_fnb option
*)
      | Diag
      | Spgmr of spils_params * 'a spils_callbacks
      | Spbcg of spils_params * 'a spils_callbacks
      | Sptfqmr of spils_params * 'a spils_callbacks
(*
      (* TODO: not for nvectors *)
      | BandedSpgmr of spils_params * bandrange
      (* TODO: not for nvectors *)
      | BandedSpbcg of spils_params * bandrange
      (* TODO: not for nvectors *)
      | BandedSptfqmr of spils_params * bandrange
*)

(*
    (* TODO: not for nvectors *)
    and 'a bdense_jac_fn = 'a bdense_jac_fn
    and 'a bband_jac_fn = 'a bband_jac_fn
    and bandrange = { mupper : int; mlower : int }
*)

    and spils_params = {
          maxl : int option;
          prec_type : Spils.preconditioning_type;
      }

    and 'a spils_callbacks =
      {
        prec_solve_fn : (('a single_tmp, 'a) jacobian_arg -> 'a prec_solve_arg
                         -> 'a -> unit) option;

        prec_setup_fn : (('a triple_tmp, 'a) jacobian_arg -> bool -> float
                         -> bool) option;

        jac_times_vec_fn : (('a single_tmp, 'a) jacobian_arg -> 'a -> 'a
                            -> unit) option;
      }

    and 'a prec_solve_arg = 'a bprec_solve_arg =
      {
        rvecB   : 'a;
        gammaB : float;
        deltaB : float;
      }

    type 'a bsession = Bsession of 'a session
    let tosession = function Bsession s -> s

    type 'a tolerance =
      | SSTolerances of float * float
      | SVTolerances of float * 'a nvector

    external ss_tolerances
        : 'a bsession -> float -> float -> unit
        = "c_cvodes_adj_ss_tolerancesb"

    external sv_tolerances
        : 'a bsession -> float -> 'a nvector -> unit
        = "c_nvec_cvodes_adj_sv_tolerancesb"

    let set_tolerances bs tol =
      match tol with
      | SSTolerances (rel, abs) -> ss_tolerances bs rel abs
      | SVTolerances (rel, abs) -> sv_tolerances bs rel abs

    external c_diag : 'a bsession -> unit
      = "c_cvode_adj_diag"

    external c_spils_set_preconditioner
      : 'a bsession -> bool -> bool -> unit
      = "c_nvec_cvode_adj_spils_set_preconditioner"

    external c_spils_spgmr
      : 'a bsession -> int -> Spils.preconditioning_type -> unit
      = "c_cvode_adj_spils_spgmr"

    external c_spils_spbcg
      : 'a bsession -> int -> Spils.preconditioning_type -> unit
      = "c_cvode_adj_spils_spbcg"

    external c_spils_sptfqmr
      : 'a bsession -> int -> Spils.preconditioning_type -> unit
      = "c_cvode_adj_spils_sptfqmr"

(* TODO serial only
    external c_spils_banded_spgmr
      : 'a session -> int -> int -> int -> int
                   -> Spils.preconditioning_type -> unit
      = "c_cvode_adj_spils_banded_spgmr"

    external c_spils_banded_spbcg
      : 'a session -> int -> int -> int -> int
                   -> Spils.preconditioning_type -> unit
      = "c_cvode_adj_spils_banded_spbcg"

    external c_spils_banded_sptfqmr
      : 'a session -> int -> int -> int -> int
                   -> Spils.preconditioning_type -> unit
      = "c_cvode_adj_spils_banded_sptfqmr"
*)

    external c_set_functional : 'a bsession -> unit
      = "c_cvode_adj_set_functional"

    let bwdsensext = function (Bsession bs) ->
      match bs.sensext with
      | BwdSensExt se -> se
      | _ -> raise AdjointNotInitialized

    (* TODO: rework for serial nvectors *)
    let set_iter_type bs iter =
      let se = bwdsensext bs in
      let optionally f = function
        | None -> ()
        | Some x -> f x
      in
      let set_precond prec_type cb =
        match prec_type with
        | Spils.PrecNone -> ()
        | Spils.PrecLeft | Spils.PrecRight | Spils.PrecBoth ->
          match cb.prec_solve_fn with
          | None -> invalid_arg "preconditioning type is not PrecNone, but no \
                                 solve function given"
          | Some solve_fn ->
            c_spils_set_preconditioner bs
              (cb.prec_setup_fn <> None)
              (cb.jac_times_vec_fn <> None);
            se.bpresolvefn <- solve_fn;
            optionally (fun f -> se.bpresetupfn <- f) cb.prec_setup_fn;
            optionally (fun f -> se.bjactimesfn <- f) cb.jac_times_vec_fn
      in
      (* Release references to all linear solver-related callbacks.  *)
      se.bpresetupfn <- dummy_bprec_setup;
      se.bpresolvefn <- dummy_bprec_solve;
      se.bjactimesfn <- dummy_bjac_times_vec;
      match iter with
      | Functional -> c_set_functional bs
      | Newton linsolv ->
        (* Iter type will be set to CV_NEWTON in the functions that set the linear
           solver.  *)
        match linsolv with
        | Diag -> c_diag bs
        | Spgmr (par, cb) ->
            let maxl = match par.maxl with None -> 0 | Some ml -> ml in
            c_spils_spgmr bs maxl par.prec_type;
            set_precond par.prec_type cb
        | Spbcg (par, cb) ->
            let maxl = match par.maxl with None -> 0 | Some ml -> ml in
            c_spils_spbcg bs maxl par.prec_type;
            set_precond par.prec_type cb
        | Sptfqmr (par, cb) ->
            let maxl = match par.maxl with None -> 0 | Some ml -> ml in
            c_spils_sptfqmr bs maxl par.prec_type;
            set_precond par.prec_type cb

    external bsession_finalize : 'a session -> unit
        = "c_cvode_bsession_finalize"

    external c_init_backward
        : 'a session Weak.t -> Cvode.lmm -> 'a iter -> 'a nvector -> float
          -> (cvode_mem * int * c_weak_ref * cvode_file)
        = "c_nvec_cvode_adj_init_backward"
    (* TODO: cvode_memB = CVodeGetAdjCVodeBmem(cvode_mem, which) *)

    external c_init_backward_1
        : 'a session Weak.t -> Cvode.lmm -> 'a iter -> 'a nvector -> float
          -> (cvode_mem * int * c_weak_ref * cvode_file)
        = "c_nvec_cvode_adj_init_backward"
    (* TODO: cvode_memB = CVodeGetAdjCVodeBmem(cvode_mem, which) *)

    let init_backward s lmm iter tol mf t0 y0 =
      let weakref = Weak.create 1 in
      let cvode_mem, which, backref, err_file =
        match mf with
        | BackBasic _ -> c_init_backward weakref lmm iter y0 t0
        | BackWithSens _ -> c_init_backward_1 weakref lmm iter y0 t0
      in
      (* cvode_mem and backref have to be immediately captured in a session and
         associated with the finalizer before we do anything else.  *)
      let bs = Bsession {
              cvode      = cvode_mem;
              backref    = backref;
              nroots     = 0;
              err_file   = err_file;

              exn_temp   = None;

              rhsfn      = (fun _ _ _ -> ());
              rootsfn    = (fun _ _ _ -> ());
              errh       = (fun _ -> ());
              errw       = (fun _ _ -> ());
              presetupfn = dummy_prec_setup;
              presolvefn = dummy_prec_solve;
              jactimesfn = dummy_jac_times_vec;

              sensext    = BwdSensExt {
                parent   = s;
                which    = which;

                brhsfn      = (match mf with
                               | BackBasic f -> f
                               | _ -> (fun _ _ _ _ -> ()));

                brhsfn1     = (match mf with
                               | BackWithSens f -> f
                               | _ -> (fun _ _ _ _ _ -> ()));

                bquadrhsfn  = (fun _ _ _ _ -> ());
                bquadrhsfn1 = (fun _ _ _ _ _ -> ());

                bpresetupfn = dummy_bprec_setup;
                bpresolvefn = dummy_bprec_solve;
                bjactimesfn = dummy_bjac_times_vec;

                (* TODO: for serial
                bdense_jac = dummy_bdense_jac;
                bband_jac = dummy_bband_jac;
                *)
              };
            } in
      Gc.finalise bsession_finalize (tosession bs);
      Weak.set weakref 0 (Some (tosession bs));
      (* Now the session is safe to use.  If any of the following fails and raises
         an exception, the GC will take care of freeing cvode_mem and backref.  *)
      set_iter_type bs iter;
      set_tolerances bs tol;
      bs

    external reinit : 'a bsession -> float -> 'a nvector -> unit
        = "c_nvec_cvodes_adj_reinit"

    external backward_normal : 'a session -> float -> unit
        = "c_cvodes_adj_backward_normal"

    external backward_one_step : 'a session -> float -> unit
        = "c_cvodes_adj_backward_one_step"

    external get : 'a bsession -> 'a nvector -> float
        = "c_nvec_cvodes_adj_backward_one_step"

    let get_dky bs = Cvode_nvector.get_dky (tosession bs)

    external set_no_sensitivity : 'a session -> unit
        = "c_cvodes_adj_set_no_sensitivity"

    external set_max_ord : 'a bsession -> int -> unit
        = "c_cvodes_set_max_ord"

    external set_max_num_steps : 'a bsession -> int -> unit
        = "c_cvodes_set_max_num_steps"

    external set_init_step : 'a bsession -> float -> unit
        = "c_cvodes_set_init_step"

    external set_min_step : 'a bsession -> float -> unit
        = "c_cvodes_set_min_step"

    external set_max_step : 'a bsession -> float -> unit
        = "c_cvodes_set_max_step"

    external set_stab_lim_det : 'a bsession -> bool -> unit
        = "c_cvodes_set_stab_lim_det"

    module Diag =
      struct
        let get_work_space bs =
          Cvode_nvector.Diag.get_work_space (tosession bs)

        let get_num_rhs_evals bs =
          Cvode_nvector.Diag.get_num_rhs_evals (tosession bs)
      end

    module Spils =
      struct
        external set_prec_type
            : 'a bsession -> Spils.preconditioning_type -> unit
            = "c_cvodes_set_prec_type"

        external set_gs_type : 'a bsession -> Spils.gramschmidt_type -> unit
            = "c_cvodes_set_gs_type"

        external set_eps_lin : 'a bsession -> float -> unit
            = "c_cvodes_set_eps_lin"

        external c_set_maxl : 'a bsession -> int -> unit
            = "c_cvodes_set_maxl"
        let set_maxl bs omaxl =
          c_set_maxl bs (match omaxl with None -> 0 | Some x -> x)

        let get_work_space bs =
          Cvode_nvector.Spils.get_work_space (tosession bs)

        let get_num_lin_iters bs =
          Cvode_nvector.Spils.get_num_lin_iters (tosession bs)

        let get_num_conv_fails bs =
          Cvode_nvector.Spils.get_num_conv_fails (tosession bs)

        let get_num_prec_evals bs =
          Cvode_nvector.Spils.get_num_prec_evals (tosession bs)

        let get_num_prec_solves bs =
          Cvode_nvector.Spils.get_num_prec_solves (tosession bs)

        let get_num_jtimes_evals bs =
          Cvode_nvector.Spils.get_num_jtimes_evals (tosession bs)

        let get_num_rhs_evals bs =
          Cvode_nvector.Spils.get_num_rhs_evals (tosession bs)
      end

    module BandPrec =
      struct
        let get_work_space bs =
          Cvode_nvector.BandPrec.get_work_space (tosession bs)
        let get_num_rhs_evals bs =
          Cvode_nvector.BandPrec.get_num_rhs_evals (tosession bs)
      end

    let get_work_space bs = Cvode_nvector.get_work_space (tosession bs)

    let get_num_steps bs = Cvode_nvector.get_num_steps (tosession bs)

    let get_num_rhs_evals bs = Cvode_nvector.get_num_rhs_evals (tosession bs)

    let get_num_lin_solv_setups bs =
      Cvode_nvector.get_num_lin_solv_setups (tosession bs)

    let get_num_err_test_fails bs =
      Cvode_nvector.get_num_err_test_fails (tosession bs)

    let get_last_order bs = Cvode_nvector.get_last_order (tosession bs)

    let get_current_order bs = Cvode_nvector.get_current_order (tosession bs)

    let get_last_step bs = Cvode_nvector.get_last_step (tosession bs)

    let get_current_step bs = Cvode_nvector.get_current_step (tosession bs)

    let get_actual_init_step bs =
      Cvode_nvector.get_actual_init_step (tosession bs)

    let get_current_time bs = Cvode_nvector.get_current_time (tosession bs)

    let get_num_stab_lim_order_reds bs =
      Cvode_nvector.get_num_stab_lim_order_reds (tosession bs)

    let get_tol_scale_factor bs =
      Cvode_nvector.get_tol_scale_factor (tosession bs)

    let get_err_weights bs = Cvode_nvector.get_err_weights (tosession bs)
    let get_est_local_errors bs =
      Cvode_nvector.get_est_local_errors (tosession bs)

    let get_integrator_stats bs =
      Cvode_nvector.get_integrator_stats (tosession bs)

    let print_integrator_stats bs =
      Cvode_nvector.print_integrator_stats (tosession bs)

    let get_num_nonlin_solv_iters bs =
      Cvode_nvector.get_num_nonlin_solv_iters (tosession bs)

    let get_num_nonlin_solv_conv_fails bs =
      Cvode_nvector.get_num_nonlin_solv_conv_fails (tosession bs)

    module Quadrature =
      struct
        type 'a _bquadrhsfn = 'a bquadrhsfn =
            QuadBasic of (float -> 'a -> 'a -> 'a -> unit)
          | QuadWithSens of (float -> 'a -> 'a array -> 'a -> 'a -> unit)
        type 'a bquadrhsfn = 'a _bquadrhsfn =
            QuadBasic of (float -> 'a -> 'a -> 'a -> unit)
          | QuadWithSens of (float -> 'a -> 'a array -> 'a -> 'a -> unit)

        external c_init : 'a bsession -> 'a nvector -> unit
            = "c_nvec_cvodes_quad_initb"
        external c_init_s : 'a bsession -> 'a nvector -> unit
            = "c_nvec_cvodes_quad_initbs"

        let init bs mf y0 =
          let se = bwdsensext bs in
          match mf with
           | QuadBasic f -> (se.bquadrhsfn <- f; c_init bs y0)
           | QuadWithSens f -> (se.bquadrhsfn1 <- f; c_init_s bs y0)

        external reinit : 'a bsession -> 'a nvector -> unit
            = "c_nvec_cvodes_quad_reinitb"

        external get : 'a bsession -> 'a nvector array -> float
            = "c_nvec_cvodes_quad_getb"

        type 'a tolerance =
            NoStepSizeControl
          | SSTolerances of float * float
          | SVTolerances of float * 'a nvector

        external set_err_con : 'a bsession -> bool -> unit
            = "c_nvec_cvodes_quad_set_err_conb"

        external sv_tolerances
            : 'a bsession -> float -> 'a nvector -> unit
            = "c_nvec_cvodes_quad_sv_tolerancesb"

        external ss_tolerances  : 'a bsession -> float -> float -> unit
            = "c_cvodes_quad_ss_tolerancesb"

        let set_tolerances bs tol =
          match tol with
          | NoStepSizeControl -> set_err_con bs false
          | SSTolerances (rel, abs) -> (set_err_con bs true;
                                        ss_tolerances bs rel abs)
          | SVTolerances (rel, abs) -> (set_err_con bs true;
                                        sv_tolerances bs rel abs)

        let get_num_rhs_evals bs =
          Quadrature.get_num_rhs_evals (tosession bs)

        let get_num_err_test_fails bs =
          Quadrature.get_num_err_test_fails (tosession bs)

        let get_quad_err_weights bs =
          Quadrature.get_quad_err_weights (tosession bs)

        let get_stats bs = Quadrature.get_stats (tosession bs)
      end
  end

