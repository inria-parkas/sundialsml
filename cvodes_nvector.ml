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

(* TODO: should we pass 'a nvectors to the get functions or just 'a's? *)

include Cvode_session_nvector

external c_alloc_nvector_array : int -> 'a array
    = "c_nvec_cvodes_alloc_nvector_array"

let add_fwdsensext s =
  match s.sensext with
  | FwdSensExt se -> ()
  | BwdSensExt _ -> failwith "Quadrature.add_fwdsensext: internal error"
  | NoSensExt ->
      s.sensext <- FwdSensExt {
        num_sensitivities = 0;
        sensarray1      = c_alloc_nvector_array 0;
        sensarray2      = c_alloc_nvector_array 0;
        quadrhsfn       = (fun _ _ _ -> ());
        senspvals       = None;
        sensrhsfn       = (fun _ _ _ _ _ _ _ -> ());
        sensrhsfn1      = (fun _ _ _ _ _ _ _ _ -> ());
        quadsensrhsfn   = (fun _ _ _ _ _ _ _ -> ());
        bsessions       = [];
      }

let num_sensitivities s =
  match s.sensext with
  | FwdSensExt se -> se.num_sensitivities
  | BwdSensExt se -> se.bnum_sensitivities
  | _ -> 0

let read_weak_fwd_ref x =
  match Weak.get x 0 with
  | Some y -> (match y.sensext with
               | FwdSensExt se -> (y, se)
               | _ -> raise (Failure "Internal error: not forward extension"))
  | None -> raise (Failure "Internal error: weak reference is dead")

let read_weak_bwd_ref x =
  match Weak.get x 0 with
  | Some y -> (match y.sensext with
               | BwdSensExt se -> (y, se)
               | _ -> raise (Failure "Internal error: not backward extension"))
  | None -> raise (Failure "Internal error: weak reference is dead")

let adjust_retcode = fun session f x ->
  try f x; 0
  with
  | Sundials.RecoverableFailure -> 1
  | e -> (session.exn_temp <- Some e; -1)

let call_quadrhsfn session t y yqdot =
  let (session, fwdsensext) = read_weak_fwd_ref session in
  adjust_retcode session (fwdsensext.quadrhsfn t y) yqdot

(* fwdsensext.sensrhsfn is called directly from C *)

let call_sensrhsfn1 session t y ydot iS yS ySdot tmp1 tmp2 =
  let (session, fwdsensext) = read_weak_fwd_ref session in
  adjust_retcode session
    (fwdsensext.sensrhsfn1 t y ydot iS yS ySdot tmp1) tmp2

(* fwdsensext.quadsensrhsfn is called directly from C *)

let call_brhsfn session t y yb ybdot =
  let (session, bwdsensext) = read_weak_bwd_ref session in
  adjust_retcode session (bwdsensext.brhsfn t y yb) ybdot

(* bwdsensext.brhsfn1 is called directly from C *)

let call_bquadrhsfn session t y yb qbdot =
  let (session, bwdsensext) = read_weak_bwd_ref session in
  adjust_retcode session (bwdsensext.bquadrhsfn t y yb) qbdot

(* bwdsensext.bquadrhsfn1 is called directly from C *)

(* the bpresetupfn is called directly from C. *)

let call_bpresolvefn session jac ps rvecB =
  let (session, bwdsensext) = read_weak_bwd_ref session in
  adjust_retcode session (bwdsensext.bpresolvefn jac ps) rvecB

let call_bjactimesfn session jac vB jvB =
  let (session, bwdsensext) = read_weak_bwd_ref session in
  adjust_retcode session (bwdsensext.bjactimesfn jac vB) jvB

let _ =
  Callback.register "c_nvec_cvodes_call_quadrhsfn"     call_quadrhsfn;
  Callback.register "c_nvec_cvodes_call_sensrhsfn1"    call_sensrhsfn1;

  Callback.register "c_nvec_cvodes_call_brhsfn"        call_brhsfn;
  Callback.register "c_nvec_cvodes_call_bquadrhsfn"    call_bquadrhsfn;
  Callback.register "c_nvec_cvodes_call_bpresolvefn"   call_bpresolvefn;
  Callback.register "c_nvec_cvodes_call_bjactimesfn"   call_bjactimesfn

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
        ("cvodes_UnrecoverableQuadRhsFuncErr",  UnrecoverableQuadRhsFuncErr);
      ]

    let fwdsensext s =
      match s.sensext with
      | FwdSensExt se -> se
      | _ -> raise QuadNotInitialized

    type 'a quadrhsfn = float -> 'a -> 'a -> unit

    external c_quad_init : 'a session -> 'a nvector -> unit
        = "c_nvec_cvodes_quad_init"

    let init s f v0 =
      add_fwdsensext s;
      let se = fwdsensext s in
      se.quadrhsfn <- f;
      c_quad_init s v0

    external reinit : 'a session -> 'a nvector -> unit
      = "c_nvec_cvodes_quad_reinit"

    external set_err_con    : 'a session -> bool -> unit
        = "c_cvodes_quad_set_err_con"
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
        = "c_nvec_cvodes_quad_get"

    external get_dky : 'a session -> float -> int -> 'a nvector -> unit
        = "c_nvec_cvodes_quad_get_dky"

    external get_num_rhs_evals       : 'a session -> int
        = "c_cvodes_quad_get_num_rhs_evals"

    external get_num_err_test_fails  : 'a session -> int
        = "c_cvodes_quad_get_num_err_test_fails"

    external get_err_weights : 'a session -> 'a nvector -> unit
        = "c_nvec_cvodes_quad_get_err_weights"

    external get_stats : 'a session -> int * int
        = "c_cvodes_quad_get_stats"
  end

module Sensitivity =
  struct
    type 'a tolerance =
        SSTolerances of float * Sundials.real_array
      | SVTolerances of float * 'a nvector array
      | EETolerances

    external set_err_con : 'a session -> bool -> unit
        = "c_cvodes_sens_set_err_con"

    external ss_tolerances  : 'a session -> float -> Sundials.real_array -> unit
        = "c_cvodes_sens_ss_tolerances"

    external ee_tolerances  : 'a session -> unit
        = "c_cvodes_sens_ee_tolerances"

    external sv_tolerances  : 'a session -> float -> 'a nvector array -> unit
        = "c_nvec_cvodes_sens_sv_tolerances"

    let set_tolerances s tol =
      let ns = num_sensitivities s in
      match tol with
      | SSTolerances (rel, abs) -> begin
            if Bigarray.Array1.dim abs <> ns
            then invalid_arg "set_tolerances: abstol has the wrong length";
            ss_tolerances s rel abs
          end
      | SVTolerances (rel, abs) -> begin
            if Array.length abs <> ns
            then invalid_arg "set_tolerances: abstol has the wrong length";
            sv_tolerances s rel abs
          end
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
        ("cvodes_BadIS",                        BadIS);
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
          (float -> 'a -> 'a -> 'a array -> 'a array -> 'a -> 'a -> unit) option
      | OneByOne of
          (float -> 'a -> 'a -> int -> 'a -> 'a -> 'a -> 'a -> unit) option
    type 'a sensrhsfn = 'a _sensrhsfn =
        AllAtOnce of
          (float -> 'a -> 'a -> 'a array -> 'a array -> 'a -> 'a -> unit) option
      | OneByOne of
          (float -> 'a -> 'a -> int -> 'a -> 'a -> 'a -> 'a -> unit) option

    type sens_params = {
        pvals  : Sundials.real_array option;
        pbar   : Sundials.real_array option;
        plist  : int array option;
      }

    let no_sens_params = { pvals = None; pbar = None; plist = None }

    external c_sens_init
        : 'a session -> sens_method -> bool -> 'a nvector array -> unit
        = "c_nvec_cvodes_sens_init"

    external c_sens_init_1
        : 'a session -> sens_method -> bool -> 'a nvector array -> unit
        = "c_nvec_cvodes_sens_init_1"

    external c_set_params : 'a session -> sens_params -> unit
        = "c_cvodes_sens_set_params"

    let set_params s ({pvals; pbar; plist} as ps) =
      let ns = num_sensitivities s in
      let np = match pvals with None -> 0 | Some p -> Bigarray.Array1.dim p in
      let check_pi v =
        if v < 0 || v >= np
        then invalid_arg "set_params: plist has an invalid entry" in
      (match pbar with
       | None -> ()
       | Some p -> if Bigarray.Array1.dim p <> ns
                   then invalid_arg "set_params: pbar has the wrong length");
      (match plist with
       | None -> ()
       | Some p -> if Array.length p <> ns
                   then invalid_arg "set_params: plist has the wrong length"
                   else Array.iter check_pi p);
      c_set_params s ps

    let init s tol fmethod sparams fm v0 =
      add_fwdsensext s;
      let se = fwdsensext s in
      let ns = Array.length v0 in
      (match fm with
       | AllAtOnce fo -> begin
           if fmethod = Staggered1 then
             failwith "Forward.init: Cannot combine AllAtOnce and Staggered1";
           (match fo with Some f -> se.sensrhsfn <- f | None -> ());
           c_sens_init s fmethod (fo <> None) v0
         end
       | OneByOne fo -> begin
           (match fo with Some f -> se.sensrhsfn1 <- f | None -> ());
           c_sens_init_1 s fmethod (fo <> None) v0
         end);
      se.num_sensitivities <- ns;
      se.senspvals <- sparams.pvals;
      se.sensarray1 <- c_alloc_nvector_array ns;
      se.sensarray2 <- c_alloc_nvector_array ns;
      set_params s sparams;
      set_tolerances s tol

    external c_reinit : 'a session -> sens_method -> 'a nvector array -> unit
        = "c_nvec_cvodes_sens_reinit"

    let reinit s sm s0 =
      if Array.length s0 <> num_sensitivities s
      then invalid_arg "reinit: wrong number of sensitivity vectors";
      c_reinit s sm s0

    external toggle_off : 'a session -> unit
        = "c_cvodes_sens_toggle_off"

    external c_get : 'a session -> 'a nvector array -> float
        = "c_nvec_cvodes_sens_get"

    let get s ys =
      if Array.length ys <> num_sensitivities s
      then invalid_arg "get: wrong number of sensitivity vectors";
      c_get s ys

    external c_get_dky : 'a session -> float -> int -> 'a nvector array -> unit
        = "c_nvec_cvodes_sens_get_dky"

    let get_dky s t k dkys =
      if Array.length dkys <> num_sensitivities s
      then invalid_arg "get_dky: wrong number of sensitivity vectors";
      c_get_dky s t k dkys

    external get1 : 'a session -> int -> 'a nvector -> float
        = "c_nvec_cvodes_sens_get1"

    external get_dky1 : 'a session -> float -> int -> int -> 'a nvector -> unit
        = "c_nvec_cvodes_sens_get_dky1"

    type dq_method = DQCentered | DQForward

    external set_dq_method : 'a session -> dq_method -> float -> unit
        = "c_cvodes_sens_set_dq_method"

    external set_max_nonlin_iters : 'a session -> int -> unit
        = "c_cvodes_sens_set_max_nonlin_iters"

    external get_num_rhs_evals : 'a session -> int
        = "c_cvodes_sens_get_num_rhs_evals"

    external get_num_rhs_evals_sens : 'a session -> int
        = "c_cvodes_sens_get_num_rhs_evals_sens"

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

    external c_get_err_weights : 'a session -> 'a nvector array -> unit
        = "c_nvec_cvodes_sens_get_err_weights"

    let get_err_weights s esweight =
      if Array.length esweight <> num_sensitivities s
      then invalid_arg "get_err_weights: wrong number of vectors";
      c_get_err_weights s esweight

    external get_num_nonlin_solv_iters : 'a session -> int
        = "c_cvodes_sens_get_num_nonlin_solv_iters"

    external get_num_nonlin_solv_conv_fails : 'a session -> int
        = "c_cvodes_sens_get_num_nonlin_solv_conv_fails"

    type nonlin_stats = {
        num_nonlin_solv_iters : int;
        num_nonlin_solv_conv_fails : int;
      }

    external get_nonlin_solv_stats : 'a session -> nonlin_stats
        = "c_cvodes_sens_get_nonlin_solv_stats"

    external c_get_num_stgr_nonlin_solv_iters
        : 'a session -> Sundials.lint_array -> unit
        = "c_cvodes_sens_get_num_stgr_nonlin_solv_iters"

    let get_num_stgr_nonlin_solv_iters s r =
      if Bigarray.Array1.dim r <> num_sensitivities s then invalid_arg
        "get_num_stgr_nonlin_solv_iters: wrong number of sensitivity vectors";
      c_get_num_stgr_nonlin_solv_iters s r

    external c_get_num_stgr_nonlin_solv_conv_fails
        : 'a session -> Sundials.lint_array -> unit
        = "c_cvodes_sens_get_num_stgr_nonlin_solv_conv_fails"

    let get_num_stgr_nonlin_solv_conv_fails s r =
      if Bigarray.Array1.dim r <> num_sensitivities s then invalid_arg
     "get_num_stgr_nonlin_solv_conv_fails: wrong number of sensitivity vectors";
      c_get_num_stgr_nonlin_solv_conv_fails s r

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
                                             UnrecoverableQuadSensRhsFuncErr);
          ]

        type 'a quadsensrhsfn =
           float          (* t *)
           -> 'a          (* y *)
           -> 'a array    (* yS *)
           -> 'a          (* yQdot *)
           -> 'a array    (* rhsvalQs *)
           -> 'a          (* tmp1 *)
           -> 'a          (* tmp2 *)
           -> unit

        external c_quadsens_init : 'a session -> 'a nvector array -> unit
            = "c_nvec_cvodes_quadsens_init"

        let init s f v0 =
          let se = fwdsensext s in
          let ns = num_sensitivities s in
          if Array.length v0 <> ns
          then invalid_arg "init: wrong number of vectors";
          se.quadsensrhsfn <- f;
          c_quadsens_init s v0

        external c_reinit : 'a session -> 'a nvector array -> unit
            = "c_nvec_cvodes_quadsens_reinit"

        let reinit s v =
          let ns = num_sensitivities s in
          if Array.length v <> ns
          then invalid_arg "reinit: wrong number of vectors";
          c_reinit s v

        type 'a tolerance =
            NoStepSizeControl
          | SSTolerances of float * Sundials.real_array
          | SVTolerances of float * 'a nvector array
          | EETolerances

        external set_err_con : 'a session -> bool -> unit
            = "c_cvodes_quadsens_set_err_con"

        external ss_tolerances  : 'a session -> float -> Sundials.real_array
                                             -> unit
            = "c_cvodes_quadsens_ss_tolerances"

        external sv_tolerances  : 'a session -> float -> 'a nvector array -> unit
            = "c_nvec_cvodes_quadsens_sv_tolerances"

        external ee_tolerances  : 'a session -> unit
            = "c_cvodes_quadsens_ee_tolerances"

        let set_tolerances s tol =
          let ns = num_sensitivities s in
          match tol with
          | NoStepSizeControl -> set_err_con s false
          | SSTolerances (rel, abs) -> begin
                if Bigarray.Array1.dim abs <> ns
                then invalid_arg "set_tolerances: abstol has the wrong length";
                set_err_con s true;
                ss_tolerances s rel abs
              end
          | SVTolerances (rel, abs) -> begin
                if Array.length abs <> ns
                then invalid_arg "set_tolerances: abstol has the wrong length";
                set_err_con s true;
                sv_tolerances s rel abs
              end
          | EETolerances -> (set_err_con s true;
                             ee_tolerances s)
    
        external c_get : 'a session -> 'a nvector array -> float
            = "c_nvec_cvodes_quadsens_get"

        let get s ys =
          let ns = num_sensitivities s in
          if Array.length ys <> ns
          then invalid_arg "get: wrong number of vectors";
          c_get s ys

        external get1 : 'a session -> int -> 'a nvector -> float
            = "c_nvec_cvodes_quadsens_get1"

        external c_get_dky : 'a session -> float -> int -> 'a nvector array
                                        -> unit
            = "c_nvec_cvodes_quadsens_get_dky"

        let get_dky s t k ys =
          let ns = num_sensitivities s in
          if Array.length ys <> ns
          then invalid_arg "get_dky: wrong number of vectors";
          c_get_dky s t k ys

        external get_dky1 : 'a session -> float -> int -> int -> 'a nvector
                                       -> unit
            = "c_nvec_cvodes_quadsens_get_dky1"

        external get_num_rhs_evals       : 'a session -> int
            = "c_cvodes_quadsens_get_num_rhs_evals"

        external get_num_err_test_fails  : 'a session -> int
            = "c_cvodes_quadsens_get_num_err_test_fails"

        external c_get_err_weights : 'a session -> 'a nvector array -> unit
            = "c_nvec_cvodes_quadsens_get_err_weights"

        let get_err_weights s esweight =
          let ns = num_sensitivities s in
          if Array.length esweight <> ns
          then invalid_arg "get_err_weights: wrong number of vectors";
          c_get_err_weights s esweight

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
    exception BadFinalTime
    exception BadOutputTime

    let _ = List.iter (fun (nm, ex) -> Callback.register_exception nm ex)
      [
        ("cvodes_AdjointNotInitialized",         AdjointNotInitialized);
        ("cvodes_NoForwardCall",                 NoForwardCall);
        ("cvodes_ForwardReinitializationFailed", ForwardReinitializationFailed);
        ("cvodes_ForwardFailed",                 ForwardFailed);
        ("cvodes_NoBackwardProblem",             NoBackwardProblem);
        ("cvodes_BadFinalTime",                  BadFinalTime);
        ("cvodes_BadOutputTime",                 BadOutputTime);
      ]

    type interpolation = IPolynomial | IHermite

    external c_init : 'a session -> int -> interpolation -> unit
        = "c_cvodes_adj_init"

    let init s nd interptype =
      add_fwdsensext s;
      c_init s nd interptype

    let fwdsensext s =
      match s.sensext with
      | FwdSensExt se -> se
      | _ -> raise AdjointNotInitialized

    external forward_normal : 'a session -> float -> 'a nvector
                                         -> float * int * Cvode.solver_result
        = "c_nvec_cvodes_adj_forward_normal"

    external forward_one_step : 'a session -> float -> 'a nvector
                                           -> float * int * Cvode.solver_result
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
        jac_u   : 'a;
        jac_ub  : 'a;
        jac_fub : 'a;
        jac_tmp : 't
      }

    type 'a iter =
      | Newton of 'a linear_solver
      | Functional

    and 'a linear_solver =
      | Diag
      | Spgmr of spils_params * 'a spils_callbacks
      | Spbcg of spils_params * 'a spils_callbacks
      | Sptfqmr of spils_params * 'a spils_callbacks

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
        rvec   : 'a;
        gamma  : float;
        delta  : float;
        left   : bool;
      }

    let spils_no_precond = {
      prec_solve_fn = None;
      prec_setup_fn = None;
      jac_times_vec_fn = None;
    }

    type 'a bsession = Bsession of 'a session
    let tosession = function Bsession s -> s

    let parent_and_which s =
      match (tosession s).sensext with
      | BwdSensExt se -> (se.parent, se.which)
      | _ -> failwith "Internal error: bsession invalid"

    type 'a tolerance =
      | SSTolerances of float * float
      | SVTolerances of float * 'a nvector

    external ss_tolerances
        : 'a session -> int -> float -> float -> unit
        = "c_cvodes_adj_ss_tolerances"

    external sv_tolerances
        : 'a session -> int -> float -> 'a nvector -> unit
        = "c_nvec_cvodes_adj_sv_tolerances"

    let set_tolerances bs tol =
      let parent, which = parent_and_which bs in
      match tol with
      | SSTolerances (rel, abs) -> ss_tolerances parent which rel abs
      | SVTolerances (rel, abs) -> sv_tolerances parent which rel abs

    external c_diag : 'a session -> int -> unit
      = "c_cvodes_adj_diag"

    external c_spils_set_preconditioner
      : 'a session -> int -> bool -> bool -> unit
      = "c_nvec_cvodes_adj_spils_set_preconditioner"

    external c_spils_spgmr
      : 'a session -> int -> int -> Spils.preconditioning_type -> unit
      = "c_cvodes_adj_spils_spgmr"

    external c_spils_spbcg
      : 'a session -> int -> int -> Spils.preconditioning_type -> unit
      = "c_cvodes_adj_spils_spbcg"

    external c_spils_sptfqmr
      : 'a session -> int -> int -> Spils.preconditioning_type -> unit
      = "c_cvodes_adj_spils_sptfqmr"

    external c_set_functional : 'a bsession -> unit
      = "c_cvodes_adj_set_functional"

    let bwdsensext = function (Bsession bs) ->
      match bs.sensext with
      | BwdSensExt se -> se
      | _ -> raise AdjointNotInitialized

    let set_iter_type bs iter =
      let se = bwdsensext bs in
      let parent, which = parent_and_which bs in
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
            c_spils_set_preconditioner parent which
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
        | Diag -> c_diag parent which
        | Spgmr (par, cb) ->
            let maxl = match par.maxl with None -> 0 | Some ml -> ml in
            c_spils_spgmr parent which maxl par.prec_type;
            set_precond par.prec_type cb
        | Spbcg (par, cb) ->
            let maxl = match par.maxl with None -> 0 | Some ml -> ml in
            c_spils_spbcg parent which maxl par.prec_type;
            set_precond par.prec_type cb
        | Sptfqmr (par, cb) ->
            let maxl = match par.maxl with None -> 0 | Some ml -> ml in
            c_spils_sptfqmr parent which maxl par.prec_type;
            set_precond par.prec_type cb

    external bsession_finalize : 'a session -> unit
        = "c_cvodes_adj_bsession_finalize"

    external c_init_backward
        : 'a session -> 'a session Weak.t
          -> (Cvode.lmm * 'a iter * float * 'a nvector)
          -> bool
          -> (cvode_mem * int * c_weak_ref * cvode_file)
        = "c_nvec_cvodes_adj_init_backward"

    let init_backward s lmm iter tol mf t0 y0 =
      let { bsessions } as se = fwdsensext s in
      let ns = num_sensitivities s in
      let weakref = Weak.create 1 in
      let cvode_mem, which, backref, err_file =
        match mf with
        | BackBasic _ -> c_init_backward s weakref (lmm, iter, t0, y0) false
        | BackWithSens _ -> c_init_backward s weakref (lmm, iter, t0, y0) true
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

                bnum_sensitivities = ns;
                bsensarray = c_alloc_nvector_array ns;

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
              };
            } in
      Gc.finalise bsession_finalize (tosession bs);
      Weak.set weakref 0 (Some (tosession bs));
      (* Now the session is safe to use.  If any of the following fails and raises
         an exception, the GC will take care of freeing cvode_mem and backref.  *)
      set_iter_type bs iter;
      set_tolerances bs tol;
      se.bsessions <- (tosession bs) :: bsessions;
      bs

    external c_reinit : 'a session -> int -> float -> 'a nvector -> unit
        = "c_nvec_cvodes_adj_reinit"

    let reinit bs tb0 yb0 =
      let parent, which = parent_and_which bs in
      c_reinit parent which tb0 yb0

    external backward_normal : 'a session -> float -> unit
        = "c_cvodes_adj_backward_normal"

    external backward_one_step : 'a session -> float -> unit
        = "c_cvodes_adj_backward_one_step"

    external c_get : 'a session -> int -> 'a nvector -> float
        = "c_nvec_cvodes_adj_get"

    let get bs yb =
      let parent, which = parent_and_which bs in
      c_get parent which yb

    let get_dky bs = Cvode_nvector.get_dky (tosession bs)

    external set_no_sensitivity : 'a session -> unit
        = "c_cvodes_adj_set_no_sensitivity"

    external c_set_max_ord : 'a session -> int -> int -> unit
        = "c_cvodes_adj_set_max_ord"

    let set_max_ord bs maxordb =
      let parent, which = parent_and_which bs in
      c_set_max_ord parent which maxordb

    external c_set_max_num_steps : 'a session -> int -> int -> unit
        = "c_cvodes_adj_set_max_num_steps"

    let set_max_num_steps bs mxstepsb =
      let parent, which = parent_and_which bs in
      c_set_max_num_steps parent which mxstepsb 

    external c_set_init_step : 'a session -> int -> float -> unit
        = "c_cvodes_adj_set_init_step"

    let set_init_step bs hinb =
      let parent, which = parent_and_which bs in
      c_set_init_step parent which hinb 

    external c_set_min_step : 'a session -> int -> float -> unit
        = "c_cvodes_adj_set_min_step"

    let set_min_step bs hminb =
      let parent, which = parent_and_which bs in
      c_set_min_step parent which hminb 

    external c_set_max_step : 'a session -> int -> float -> unit
        = "c_cvodes_adj_set_max_step"

    let set_max_step bs hmaxb =
      let parent, which = parent_and_which bs in
      c_set_max_step parent which hmaxb 

    external c_set_stab_lim_det : 'a session -> int -> bool -> unit
        = "c_cvodes_adj_set_stab_lim_det"

    let set_stab_lim_det bs stldetb =
      let parent, which = parent_and_which bs in
      c_set_stab_lim_det parent which stldetb

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
            = "c_cvodes_adj_spils_set_prec_type"

        external set_gs_type : 'a bsession -> Spils.gramschmidt_type -> unit
            = "c_cvodes_adj_spils_set_gs_type"

        external set_eps_lin : 'a bsession -> float -> unit
            = "c_cvodes_adj_spils_set_eps_lin"

        external c_set_maxl : 'a bsession -> int -> unit
            = "c_cvodes_adj_spils_set_maxl"

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

        external c_quad_initb : 'a session -> int -> 'a nvector -> unit
            = "c_nvec_cvodes_adjquad_initb"
        external c_quad_initbs : 'a session -> int -> 'a nvector -> unit
            = "c_nvec_cvodes_adjquad_initbs"

        let init bs mf y0 =
          let parent, which = parent_and_which bs in
          let se = bwdsensext bs in
          match mf with
           | QuadBasic f -> (se.bquadrhsfn <- f;
                             c_quad_initb parent which y0)
           | QuadWithSens f -> (se.bquadrhsfn1 <- f;
                                c_quad_initbs parent which y0)

        external c_reinit : 'a session -> int -> 'a nvector -> unit
            = "c_nvec_cvodes_adjquad_reinit"

        let reinit bs yqb0 =
          let parent, which = parent_and_which bs in
          c_reinit parent which yqb0

        external c_get : 'a session -> int -> 'a nvector -> float
            = "c_nvec_cvodes_adjquad_get"

        let get bs yqb =
          let parent, which = parent_and_which bs in
          c_get parent which yqb

        type 'a tolerance =
            NoStepSizeControl
          | SSTolerances of float * float
          | SVTolerances of float * 'a nvector

        external set_err_con : 'a bsession -> bool -> unit
            = "c_cvodes_adjquad_set_err_con"

        external sv_tolerances
            : 'a session -> int -> float -> 'a nvector -> unit
            = "c_nvec_cvodes_adjquad_sv_tolerances"

        external ss_tolerances  : 'a session -> int -> float -> float -> unit
            = "c_cvodes_adjquad_ss_tolerances"

        let set_tolerances bs tol =
          let parent, which = parent_and_which bs in
          match tol with
          | NoStepSizeControl -> set_err_con bs false
          | SSTolerances (rel, abs) -> (set_err_con bs true;
                                        ss_tolerances parent which rel abs)
          | SVTolerances (rel, abs) -> (set_err_con bs true;
                                        sv_tolerances parent which rel abs)

        let get_num_rhs_evals bs =
          Quadrature.get_num_rhs_evals (tosession bs)

        let get_num_err_test_fails bs =
          Quadrature.get_num_err_test_fails (tosession bs)

        let get_err_weights bs =
          Quadrature.get_err_weights (tosession bs)

        let get_stats bs = Quadrature.get_stats (tosession bs)
      end
  end

