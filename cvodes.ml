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

external c_alloc_nvector_array : int -> 'a array
    = "c_cvodes_alloc_nvector_array"

let add_fwdsensext s =
  match s.sensext with
  | FwdSensExt se -> ()
  | BwdSensExt _ -> failwith "Quadrature.add_fwdsensext: internal error"
  | NoSensExt ->
      s.sensext <- FwdSensExt {
        num_sensitivities = 0;
        sensarray1      = c_alloc_nvector_array 0;
        sensarray2      = c_alloc_nvector_array 0;
        quadrhsfn       = dummy_quadrhsfn;
        senspvals       = None;
        sensrhsfn       = dummy_sensrhsfn;
        sensrhsfn1      = dummy_sensrhsfn1;
        quadsensrhsfn   = dummy_quadsensrhsfn;
        bsessions       = [];
      }

let num_sensitivities s =
  match s.sensext with
  | FwdSensExt se -> se.num_sensitivities
  | BwdSensExt se -> se.bnum_sensitivities
  | _ -> 0

module Quadrature =
  struct
    include QuadratureTypes

    exception QuadNotInitialized
    exception QuadRhsFuncFailure
    exception FirstQuadRhsFuncFailure
    exception RepeatedQuadRhsFuncFailure
    exception UnrecoverableQuadRhsFuncFailure

    let fwdsensext s =
      match s.sensext with
      | FwdSensExt se -> se
      | _ -> raise QuadNotInitialized

    external c_quad_init : ('a, 'k) session -> ('a, 'k) nvector -> unit
        = "c_cvodes_quad_init"

    let init s f v0 =
      add_fwdsensext s;
      let se = fwdsensext s in
      se.quadrhsfn <- f;
      c_quad_init s v0

    external reinit : ('a, 'k) session -> ('a, 'k) nvector -> unit
      = "c_cvodes_quad_reinit"

    external set_err_con    : ('a, 'k) session -> bool -> unit
        = "c_cvodes_quad_set_err_con"
    external sv_tolerances
        : ('a, 'k) session -> float -> ('a, 'k) nvector -> unit
        = "c_cvodes_quad_sv_tolerances"
    external ss_tolerances  : ('a, 'k) session -> float -> float -> unit
        = "c_cvodes_quad_ss_tolerances"

    type ('a, 'k) tolerance =
        NoStepSizeControl
      | SStolerances of float * float
      | SVtolerances of float * ('a, 'k) nvector

    let set_tolerances s tol =
      match tol with
      | NoStepSizeControl -> set_err_con s false
      | SStolerances (rel, abs) -> (ss_tolerances s rel abs;
                                    set_err_con s true)
      | SVtolerances (rel, abs) -> (sv_tolerances s rel abs;
                                    set_err_con s true)

    external get : ('a, 'k) session -> ('a, 'k) nvector -> float
        = "c_cvodes_quad_get"

    external get_dky
        : ('a, 'k) session -> float -> int -> ('a, 'k) nvector -> unit
        = "c_cvodes_quad_get_dky"

    external get_num_rhs_evals       : ('a, 'k) session -> int
        = "c_cvodes_quad_get_num_rhs_evals"

    external get_num_err_test_fails  : ('a, 'k) session -> int
        = "c_cvodes_quad_get_num_err_test_fails"

    external get_err_weights : ('a, 'k) session -> ('a, 'k) nvector -> unit
        = "c_cvodes_quad_get_err_weights"

    external get_stats : ('a, 'k) session -> int * int
        = "c_cvodes_quad_get_stats"
  end

module Sensitivity =
  struct
    include SensitivityTypes

    type ('a, 'k) tolerance =
        SStolerances of float * Sundials.RealArray.t
      | SVtolerances of float * ('a, 'k) nvector array
      | EEtolerances

    external set_err_con : ('a, 'k) session -> bool -> unit
        = "c_cvodes_sens_set_err_con"

    external ss_tolerances
        : ('a, 'k) session -> float -> Sundials.RealArray.t -> unit
        = "c_cvodes_sens_ss_tolerances"

    external ee_tolerances  : ('a, 'k) session -> unit
        = "c_cvodes_sens_ee_tolerances"

    external sv_tolerances
        : ('a, 'k) session -> float -> ('a, 'k) nvector array -> unit
        = "c_cvodes_sens_sv_tolerances"

    let set_tolerances s tol =
      let ns = num_sensitivities s in
      match tol with
      | SStolerances (rel, abs) -> begin
            if Bigarray.Array1.dim abs <> ns
            then invalid_arg "set_tolerances: abstol has the wrong length";
            ss_tolerances s rel abs
          end
      | SVtolerances (rel, abs) -> begin
            if Array.length abs <> ns
            then invalid_arg "set_tolerances: abstol has the wrong length";
            sv_tolerances s rel abs
          end
      | EEtolerances -> ee_tolerances s

    exception SensNotInitialized
    exception SensRhsFuncFailure
    exception FirstSensRhsFuncFailure
    exception RepeatedSensRhsFuncFailure
    exception UnrecoverableSensRhsFuncFailure
    exception BadSensIdentifier

    let fwdsensext s =
      match s.sensext with
      | FwdSensExt se -> se
      | _ -> raise SensNotInitialized

    type sens_method =
        Simultaneous
      | Staggered
      | Staggered1

    type sens_params = {
        pvals  : Sundials.RealArray.t option;
        pbar   : Sundials.RealArray.t option;
        plist  : int array option;
      }

    let no_sens_params = { pvals = None; pbar = None; plist = None }

    external c_sens_init : ('a, 'k) session -> sens_method -> bool
                                -> ('a, 'k) nvector array -> unit
        = "c_cvodes_sens_init"

    external c_sens_init_1 : ('a, 'k) session -> sens_method -> bool
                                -> ('a, 'k) nvector array -> unit
        = "c_cvodes_sens_init_1"

    external c_set_params : ('a, 'k) session -> sens_params -> unit
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

    external c_reinit
        : ('a, 'k) session -> sens_method -> ('a, 'k) nvector array -> unit
        = "c_cvodes_sens_reinit"

    let reinit s sm s0 =
      if Array.length s0 <> num_sensitivities s
      then invalid_arg "reinit: wrong number of sensitivity vectors";
      c_reinit s sm s0

    external toggle_off : ('a, 'k) session -> unit
        = "c_cvodes_sens_toggle_off"

    external c_get : ('a, 'k) session -> ('a, 'k) nvector array -> float
        = "c_cvodes_sens_get"

    let get s ys =
      if Array.length ys <> num_sensitivities s
      then invalid_arg "get: wrong number of sensitivity vectors";
      c_get s ys

    external c_get_dky
        : ('a, 'k) session -> float -> int -> ('a, 'k) nvector array -> unit
        = "c_cvodes_sens_get_dky"

    let get_dky s t k dkys =
      if Array.length dkys <> num_sensitivities s
      then invalid_arg "get_dky: wrong number of sensitivity vectors";
      c_get_dky s t k dkys

    external get1 : ('a, 'k) session -> int -> ('a, 'k) nvector -> float
        = "c_cvodes_sens_get1"

    external get_dky1
        : ('a, 'k) session -> float -> int -> int -> ('a, 'k) nvector -> unit
        = "c_cvodes_sens_get_dky1"

    type dq_method = DQCentered | DQForward

    external set_dq_method : ('a, 'k) session -> dq_method -> float -> unit
        = "c_cvodes_sens_set_dq_method"

    external set_max_nonlin_iters : ('a, 'k) session -> int -> unit
        = "c_cvodes_sens_set_max_nonlin_iters"

    external get_num_rhs_evals : ('a, 'k) session -> int
        = "c_cvodes_sens_get_num_rhs_evals"

    external get_num_rhs_evals_sens : ('a, 'k) session -> int
        = "c_cvodes_sens_get_num_rhs_evals_sens"

    external get_num_err_test_fails : ('a, 'k) session -> int
        = "c_cvodes_sens_get_num_err_test_fails"

    external get_num_lin_solv_setups : ('a, 'k) session -> int
        = "c_cvodes_sens_get_num_lin_solv_setups"

    type sensitivity_stats = {
        num_rhs_evals : int;
        num_sens_evals :int;
        num_err_test_fails : int;
        num_lin_solv_setups :int;
      }

    external get_stats : ('a, 'k) session -> sensitivity_stats
        = "c_cvodes_sens_get_stats"

    external c_get_err_weights
        : ('a, 'k) session -> ('a, 'k) nvector array -> unit
        = "c_cvodes_sens_get_err_weights"

    let get_err_weights s esweight =
      if Array.length esweight <> num_sensitivities s
      then invalid_arg "get_err_weights: wrong number of vectors";
      c_get_err_weights s esweight

    external get_num_nonlin_solv_iters : ('a, 'k) session -> int
        = "c_cvodes_sens_get_num_nonlin_solv_iters"

    external get_num_nonlin_solv_conv_fails : ('a, 'k) session -> int
        = "c_cvodes_sens_get_num_nonlin_solv_conv_fails"

    external get_nonlin_solv_stats : ('a, 'k) session -> int * int
        = "c_cvodes_sens_get_nonlin_solv_stats"

    external c_get_num_stgr_nonlin_solv_iters
        : ('a, 'k) session -> Sundials.LintArray.t -> unit
        = "c_cvodes_sens_get_num_stgr_nonlin_solv_iters"

    let get_num_stgr_nonlin_solv_iters s r =
      if Bigarray.Array1.dim r <> num_sensitivities s then invalid_arg
        "get_num_stgr_nonlin_solv_iters: wrong number of sensitivity vectors";
      c_get_num_stgr_nonlin_solv_iters s r

    external c_get_num_stgr_nonlin_solv_conv_fails
        : ('a, 'k) session -> Sundials.LintArray.t -> unit
        = "c_cvodes_sens_get_num_stgr_nonlin_solv_conv_fails"

    let get_num_stgr_nonlin_solv_conv_fails s r =
      if Bigarray.Array1.dim r <> num_sensitivities s then invalid_arg
     "get_num_stgr_nonlin_solv_conv_fails: wrong number of sensitivity vectors";
      c_get_num_stgr_nonlin_solv_conv_fails s r

    module Quadrature =
      struct
        include QuadratureTypes

        exception QuadSensNotInitialized
        exception QuadSensRhsFuncFailure
        exception FirstQuadSensRhsFuncFailure
        exception RepeatedQuadSensRhsFuncFailure
        exception UnrecoverableQuadSensRhsFuncFailure

        external c_quadsens_init
            : ('a, 'k) session -> bool -> ('a, 'k) nvector array -> unit
            = "c_cvodes_quadsens_init"

        let init s ?fQS v0 =
          let se = fwdsensext s in
          let ns = num_sensitivities s in
          if Array.length v0 <> ns
          then invalid_arg "init: wrong number of vectors";
          match fQS with
          | Some f -> se.quadsensrhsfn <- f;
                      c_quadsens_init s true v0
          | None -> c_quadsens_init s false v0

        external c_reinit : ('a, 'k) session -> ('a, 'k) nvector array -> unit
            = "c_cvodes_quadsens_reinit"

        let reinit s v =
          let ns = num_sensitivities s in
          if Array.length v <> ns
          then invalid_arg "reinit: wrong number of vectors";
          c_reinit s v

        type ('a, 'k) tolerance =
            NoStepSizeControl
          | SStolerances of float * Sundials.RealArray.t
          | SVtolerances of float * ('a, 'k) nvector array
          | EEtolerances

        external set_err_con : ('a, 'k) session -> bool -> unit
            = "c_cvodes_quadsens_set_err_con"

        external ss_tolerances
            : ('a, 'k) session -> float -> Sundials.RealArray.t -> unit
            = "c_cvodes_quadsens_ss_tolerances"

        external sv_tolerances
            : ('a, 'k) session -> float -> ('a, 'k) nvector array -> unit
            = "c_cvodes_quadsens_sv_tolerances"

        external ee_tolerances  : ('a, 'k) session -> unit
            = "c_cvodes_quadsens_ee_tolerances"

        let set_tolerances s tol =
          let ns = num_sensitivities s in
          match tol with
          | NoStepSizeControl -> set_err_con s false
          | SStolerances (rel, abs) -> begin
                if Bigarray.Array1.dim abs <> ns
                then invalid_arg "set_tolerances: abstol has the wrong length";
                ss_tolerances s rel abs;
                set_err_con s true
              end
          | SVtolerances (rel, abs) -> begin
                if Array.length abs <> ns
                then invalid_arg "set_tolerances: abstol has the wrong length";
                sv_tolerances s rel abs;
                set_err_con s true
              end
          | EEtolerances -> (ee_tolerances s;
                             set_err_con s true)
    
        external c_get : ('a, 'k) session -> ('a, 'k) nvector array -> float
            = "c_cvodes_quadsens_get"

        let get s ys =
          let ns = num_sensitivities s in
          if Array.length ys <> ns
          then invalid_arg "get: wrong number of vectors";
          c_get s ys

        external get1 : ('a, 'k) session -> int -> ('a, 'k) nvector -> float
            = "c_cvodes_quadsens_get1"

        external c_get_dky
            : ('a, 'k) session -> float -> int -> ('a, 'k) nvector array -> unit
            = "c_cvodes_quadsens_get_dky"

        let get_dky s t k ys =
          let ns = num_sensitivities s in
          if Array.length ys <> ns
          then invalid_arg "get_dky: wrong number of vectors";
          c_get_dky s t k ys

        external get_dky1 : ('a, 'k) session -> float -> int -> int
                                      -> ('a, 'k) nvector -> unit
            = "c_cvodes_quadsens_get_dky1"

        external get_num_rhs_evals       : ('a, 'k) session -> int
            = "c_cvodes_quadsens_get_num_rhs_evals"

        external get_num_err_test_fails  : ('a, 'k) session -> int
            = "c_cvodes_quadsens_get_num_err_test_fails"

        external c_get_err_weights
            : ('a, 'k) session -> ('a, 'k) nvector array -> unit
            = "c_cvodes_quadsens_get_err_weights"

        let get_err_weights s esweight =
          let ns = num_sensitivities s in
          if Array.length esweight <> ns
          then invalid_arg "get_err_weights: wrong number of vectors";
          c_get_err_weights s esweight

        external get_stats : ('a, 'k) session -> int * int
            = "c_cvodes_quadsens_get_stats"
      end
    end

module Adjoint =
  struct
    include AdjointTypes

    exception AdjointNotInitialized
    exception NoForwardCall
    exception ForwardReinitFailure
    exception ForwardFailure
    exception NoBackwardProblem
    exception BadFinalTime
    exception BadOutputTime

    let optionally f = function
      | None -> ()
      | Some x -> f x

    type interpolation = IPolynomial | IHermite

    external c_init : ('a, 'k) session -> int -> interpolation -> unit
        = "c_cvodes_adj_init"

    let init s nd interptype =
      add_fwdsensext s;
      c_init s nd interptype

    let fwdsensext s =
      match s.sensext with
      | FwdSensExt se -> se
      | _ -> raise AdjointNotInitialized

    external forward_normal : ('a, 'k) session -> float -> ('a, 'k) nvector
                                         -> float * int * Cvode.solver_result
        = "c_cvodes_adj_forward_normal"

    external forward_one_step : ('a, 'k) session -> float -> ('a, 'k) nvector
                                         -> float * int * Cvode.solver_result
        = "c_cvodes_adj_forward_one_step"

    type 'a triple = 'a * 'a * 'a

    type bandrange = Cvode_impl.bandrange = { mupper : int; mlower : int; }

    type ('data, 'kind) iter =
      | Newton of ('data, 'kind) linear_solver
      | Functional

    let parent_and_which s =
      match (tosession s).sensext with
      | BwdSensExt se -> (se.parent, se.which)
      | _ -> failwith "Internal error: bsession invalid"

    type ('a, 'k) tolerance =
      | SStolerances of float * float
      | SVtolerances of float * ('a, 'k) nvector

    external ss_tolerances
        : ('a, 'k) session -> int -> float -> float -> unit
        = "c_cvodes_adj_ss_tolerances"

    external sv_tolerances
        : ('a, 'k) session -> int -> float -> ('a, 'k) nvector -> unit
        = "c_cvodes_adj_sv_tolerances"

    let set_tolerances bs tol =
      let parent, which = parent_and_which bs in
      match tol with
      | SStolerances (rel, abs) -> ss_tolerances parent which rel abs
      | SVtolerances (rel, abs) -> sv_tolerances parent which rel abs

    external c_set_functional : ('a, 'k) session -> int -> unit
      = "c_cvodes_adj_set_functional"

    let bwdsensext = function (Bsession bs) ->
      match bs.sensext with
      | BwdSensExt se -> se
      | _ -> raise AdjointNotInitialized

    let set_iter_type bs iter nv =
      match iter with
      | Functional ->
          let parent, which = parent_and_which bs in
          c_set_functional parent which
      | Newton linsolv -> linsolv bs nv
        (* Iter type will be set to CV_NEWTON in the functions that
           set the linear solver.  *)

    external backward_normal : ('a, 'k) session -> float -> unit
        = "c_cvodes_adj_backward_normal"

    external backward_one_step : ('a, 'k) session -> float -> unit
        = "c_cvodes_adj_backward_one_step"

    external c_get : ('a, 'k) session -> int -> ('a, 'k) nvector -> float
        = "c_cvodes_adj_get"

    let get bs yb =
      let parent, which = parent_and_which bs in
      c_get parent which yb

    let get_dky bs = Cvode.get_dky (tosession bs)

    external set_no_sensitivity : ('a, 'k) session -> unit
        = "c_cvodes_adj_set_no_sensitivity"

    external c_set_max_ord : ('a, 'k) session -> int -> int -> unit
        = "c_cvodes_adj_set_max_ord"

    let set_max_ord bs maxordb =
      let parent, which = parent_and_which bs in
      c_set_max_ord parent which maxordb

    external c_set_max_num_steps : ('a, 'k) session -> int -> int -> unit
        = "c_cvodes_adj_set_max_num_steps"

    let set_max_num_steps bs mxstepsb =
      let parent, which = parent_and_which bs in
      c_set_max_num_steps parent which mxstepsb 

    external c_set_init_step : ('a, 'k) session -> int -> float -> unit
        = "c_cvodes_adj_set_init_step"

    let set_init_step bs hinb =
      let parent, which = parent_and_which bs in
      c_set_init_step parent which hinb 

    external c_set_min_step : ('a, 'k) session -> int -> float -> unit
        = "c_cvodes_adj_set_min_step"

    let set_min_step bs hminb =
      let parent, which = parent_and_which bs in
      c_set_min_step parent which hminb 

    external c_set_max_step : ('a, 'k) session -> int -> float -> unit
        = "c_cvodes_adj_set_max_step"

    let set_max_step bs hmaxb =
      let parent, which = parent_and_which bs in
      c_set_max_step parent which hmaxb 

    external c_set_stab_lim_det : ('a, 'k) session -> int -> bool -> unit
        = "c_cvodes_adj_set_stab_lim_det"

    let set_stab_lim_det bs stldetb =
      let parent, which = parent_and_which bs in
      c_set_stab_lim_det parent which stldetb

    module Diag =
      struct
        external c_diag : ('a, 'k) session -> int -> unit
          = "c_cvodes_adj_diag"

        let solver bs _ =
          let parent, which = parent_and_which bs in
          c_diag parent which

        let get_work_space bs =
          Cvode.Diag.get_work_space (tosession bs)

        let get_num_rhs_evals bs =
          Cvode.Diag.get_num_rhs_evals (tosession bs)
      end

    module Dls =
      struct
        include DlsTypes

        external c_dls_dense : serial_session -> int -> int -> bool -> unit
          = "c_cvodes_adj_dls_dense"

        external c_dls_lapack_dense
          : serial_session -> int -> int -> bool -> unit
          = "c_cvodes_adj_dls_lapack_dense"

        external c_dls_band
          : (serial_session * int) -> int -> int -> int -> bool -> unit
          = "c_cvodes_adj_dls_band"

        external c_dls_lapack_band
          : (serial_session * int) -> int -> int -> int -> bool -> unit
          = "c_cvodes_adj_dls_lapack_band"

        let dense ?jac () bs nv =
          let parent, which = parent_and_which bs in
          let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
          c_dls_dense parent which neqs (jac <> None);
          (tosession bs).ls_callbacks <-
            match jac with
            | None -> NoCallbacks
            | Some f -> BDenseCallback { jacfn = f; dmat = None }

        let lapack_dense ?jac () bs nv =
          let parent, which = parent_and_which bs in
          let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
          c_dls_lapack_dense parent which neqs (jac <> None);
          (tosession bs).ls_callbacks <-
            match jac with
            | None -> NoCallbacks
            | Some f -> BDenseCallback { jacfn = f; dmat = None }

        type ('data, 'kind) linear_solver =
          ('data, 'kind) bsession -> ('data, 'kind) nvector -> unit

        let band ?jac p bs nv =
          let parent, which = parent_and_which bs in
          let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
          c_dls_band (parent, which) neqs p.mupper p.mlower (jac <> None);
          (tosession bs).ls_callbacks <-
            match jac with
            | None -> NoCallbacks
            | Some f -> BBandCallback { bjacfn = f; bmat = None }

        let lapack_band ?jac p bs nv =
          let parent, which = parent_and_which bs in
          let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
          c_dls_lapack_band (parent,which) neqs p.mupper p.mlower (jac <> None);
          (tosession bs).ls_callbacks <-
            match jac with
            | None -> NoCallbacks
            | Some f -> BBandCallback { bjacfn = f; bmat = None }

        let invalidate_callback (type d) (type k) (s : (d, k) session) =
          match s.ls_callbacks with
          | BDenseCallback ({ dmat = Some d } as cb) ->
              Dls.DenseMatrix.invalidate d;
              cb.dmat <- None
          | BBandCallback  ({ bmat = Some d } as cb) ->
              Dls.BandMatrix.invalidate d;
              cb.bmat <- None
          | _ -> ()
      end

    module Spils =
      struct
        include SpilsTypes

        external c_set_preconditioner
          : ('a, 'k) session -> int -> bool -> unit
          = "c_cvodes_adj_spils_set_preconditioner"

        external c_set_jac_times_vec_fn
          : ('a, 'k) session -> int -> bool -> unit
          = "c_cvodes_adj_spils_set_jac_times_vec_fn"

        external c_spgmr
          : ('a, 'k) session -> int -> int -> Spils.preconditioning_type -> unit
          = "c_cvodes_adj_spils_spgmr"

        external c_spbcg
          : ('a, 'k) session -> int -> int -> Spils.preconditioning_type -> unit
          = "c_cvodes_adj_spils_spbcg"

        external c_sptfqmr
          : ('a, 'k) session -> int -> int -> Spils.preconditioning_type -> unit
          = "c_cvodes_adj_spils_sptfqmr"

        let init_preconditioner solve setup jac_times bs parent which nv =
          c_set_preconditioner parent which (setup <> None);
          c_set_jac_times_vec_fn parent which (jac_times <> None);
          (tosession bs).ls_callbacks <-
            BSpilsCallback { prec_solve_fn = solve;
                             prec_setup_fn = setup;
                             jac_times_vec_fn = jac_times }

        let prec_none = InternalPrecNone
        let prec_left ?setup ?jac_times_vec solve =
          InternalPrecLeft (init_preconditioner solve setup jac_times_vec)
        let prec_right ?setup ?jac_times_vec solve =
          InternalPrecRight (init_preconditioner solve setup jac_times_vec)
        let prec_both ?setup ?jac_times_vec solve =
          InternalPrecBoth (init_preconditioner solve setup jac_times_vec)

        let init_spils init maxl prec bs nv =
          let parent, which = parent_and_which bs in
          let with_prec prec_type set_prec =
            init parent which maxl prec_type;
            set_prec bs parent which nv
          in
          match prec with
          | InternalPrecNone -> init parent which maxl Spils.PrecNone
          | InternalPrecLeft set_prec  -> with_prec Spils.PrecLeft set_prec
          | InternalPrecRight set_prec -> with_prec Spils.PrecRight set_prec
          | InternalPrecBoth set_prec  -> with_prec Spils.PrecBoth set_prec

        let spgmr ?(maxl=0) prec bs nv =
          init_spils c_spgmr maxl prec bs nv

        let spbcg ?(maxl=0) prec bs nv =
          init_spils c_spbcg maxl prec bs nv

        let sptfqmr ?(maxl=0) prec bs nv =
          init_spils c_sptfqmr maxl prec bs nv

        external c_set_prec_type
            : ('a, 'k) bsession -> Spils.preconditioning_type -> unit
            = "c_cvodes_adj_spils_set_prec_type"

        let set_prec_type s = function
          | PrecNone -> c_set_prec_type s Spils.PrecNone
          | PrecLeft -> c_set_prec_type s Spils.PrecLeft
          | PrecRight -> c_set_prec_type s Spils.PrecRight
          | PrecBoth -> c_set_prec_type s Spils.PrecBoth

        let set_preconditioner bs ?setup solve =
          match (tosession bs).ls_callbacks with
          | BSpilsCallback cbs ->
            let parent, which = parent_and_which bs in
            c_set_preconditioner parent which (setup <> None);
            (tosession bs).ls_callbacks <-
              BSpilsCallback { cbs with
                               prec_setup_fn = setup;
                               prec_solve_fn = solve }
          | BSpilsBandCallback _ ->
            failwith "User-defined preconditioner not in use"
          | _ -> failwith "spils solver not in use"

        let set_jac_times_vec_fn bs f =
          match (tosession bs).ls_callbacks with
          | BSpilsCallback cbs ->
            let parent, which = parent_and_which bs in
            c_set_jac_times_vec_fn parent which true;
            (tosession bs).ls_callbacks <-
              BSpilsCallback { cbs with jac_times_vec_fn = Some f }
          | BSpilsBandCallback _ ->
            let parent, which = parent_and_which bs in
            c_set_jac_times_vec_fn parent which false;
            (tosession bs).ls_callbacks <- BSpilsBandCallback (Some f)
          | _ -> failwith "spils solver not in use"

        let clear_jac_times_vec_fn bs =
          match (tosession bs).ls_callbacks with
          | BSpilsCallback cbs ->
            let parent, which = parent_and_which bs in
            c_set_jac_times_vec_fn parent which false;
            (tosession bs).ls_callbacks <-
              BSpilsCallback { cbs with jac_times_vec_fn = None }
          | BSpilsBandCallback _ ->
            let parent, which = parent_and_which bs in
            c_set_jac_times_vec_fn parent which false;
            (tosession bs).ls_callbacks <- BSpilsBandCallback None
          | _ -> failwith "spils solver not in use"

        external set_gs_type
            : ('a, 'k) bsession -> gramschmidt_type -> unit
            = "c_cvodes_adj_spils_set_gs_type"

        external set_eps_lin : ('a, 'k) bsession -> float -> unit
            = "c_cvodes_adj_spils_set_eps_lin"

        external c_set_maxl : ('a, 'k) bsession -> int -> unit
            = "c_cvodes_adj_spils_set_maxl"

        let set_maxl bs omaxl =
          c_set_maxl bs (match omaxl with None -> 0 | Some x -> x)

        let get_work_space bs =
          Cvode.Spils.get_work_space (tosession bs)

        let get_num_lin_iters bs =
          Cvode.Spils.get_num_lin_iters (tosession bs)

        let get_num_conv_fails bs =
          Cvode.Spils.get_num_conv_fails (tosession bs)

        let get_num_prec_evals bs =
          Cvode.Spils.get_num_prec_evals (tosession bs)

        let get_num_prec_solves bs =
          Cvode.Spils.get_num_prec_solves (tosession bs)

        let get_num_jtimes_evals bs =
          Cvode.Spils.get_num_jtimes_evals (tosession bs)

        let get_num_rhs_evals bs =
          Cvode.Spils.get_num_rhs_evals (tosession bs)

        module Banded = struct
          external c_set_preconditioner
            : ('a, 'k) session -> int -> int -> int -> int -> unit
            = "c_cvodes_adj_spils_set_preconditioner"

          let init_preconditioner jac_times_vec bandrange bs parent which nv =
            c_set_preconditioner parent which
              (RealArray.length (Nvector.unwrap nv))
              bandrange.mupper bandrange.mlower;
            c_set_jac_times_vec_fn parent which (jac_times_vec <> None);
            (tosession bs).ls_callbacks <- BSpilsBandCallback jac_times_vec

          let prec_none = InternalPrecNone
          let prec_left ?jac_times_vec bandrange =
            InternalPrecLeft (init_preconditioner jac_times_vec bandrange)
          let prec_right ?jac_times_vec bandrange =
            InternalPrecRight (init_preconditioner jac_times_vec bandrange)
          let prec_both ?jac_times_vec bandrange =
            InternalPrecBoth (init_preconditioner jac_times_vec bandrange)

          let get_work_space bs =
            Cvode.Spils.Banded.get_work_space (tosession bs)

          let get_num_rhs_evals bs =
            Cvode.Spils.Banded.get_num_rhs_evals (tosession bs)
        end
      end

    external c_bsession_finalize : ('a, 'k) session -> unit
        = "c_cvodes_adj_bsession_finalize"

    let bsession_finalize s =
      Dls.invalidate_callback s;
      c_bsession_finalize s

    external c_init_backward
        : ('a, 'k) session -> ('a, 'k) session Weak.t
          -> (Cvode.lmm * ('a, 'k) iter * float * ('a, 'k) nvector)
          -> bool
          -> (cvode_mem * int * c_weak_ref * cvode_file)
        = "c_cvodes_adj_init_backward"

    let init_backward s lmm iter tol mf t0 y0 =
      let { bsessions } as se = fwdsensext s in
      let ns = num_sensitivities s in
      let checkvec = Nvector.check y0 in
      let weakref = Weak.create 1 in
      let cvode_mem, which, backref, err_file =
        match mf with
        | NoSens _ -> c_init_backward s weakref (lmm, iter, t0, y0) false
        | WithSens _ -> c_init_backward s weakref (lmm, iter, t0, y0) true
      in
      (* cvode_mem and backref have to be immediately captured in a session and
         associated with the finalizer before we do anything else.  *)
      let bs = Bsession {
              cvode        = cvode_mem;
              backref      = backref;
              nroots       = 0;
              err_file     = err_file;
              checkvec     = checkvec;

              exn_temp     = None;

              rhsfn        = dummy_rhsfn;
              rootsfn      = dummy_rootsfn;
              errh         = dummy_errh;
              errw         = dummy_errw;
              ls_callbacks = NoCallbacks;

              sensext    = BwdSensExt {
                parent   = s;
                which    = which;

                bnum_sensitivities = ns;
                bsensarray = c_alloc_nvector_array ns;

                brhsfn      = (match mf with
                               | NoSens f -> f
                               | _ -> dummy_brhsfn);

                brhsfn_sens = (match mf with
                               | WithSens f -> f
                               | _ -> dummy_brhsfn_sens);

                bquadrhsfn  = dummy_bquadrhsfn;
                bquadrhsfn_sens = dummy_bquadrhsfn_sens;
              };
            } in
      Gc.finalise bsession_finalize (tosession bs);
      Weak.set weakref 0 (Some (tosession bs));
      (* Now the session is safe to use.  If any of the following fails and
         raises an exception, the GC will take care of freeing cvode_mem and
         backref. *)
      set_iter_type bs iter y0;
      set_tolerances bs tol;
      se.bsessions <- (tosession bs) :: bsessions;
      bs

    external c_reinit
        : ('a, 'k) session -> int -> float -> ('a, 'k) nvector -> unit
        = "c_cvodes_adj_reinit"

    let reinit bs ?iter_type tb0 yb0 =
      let parent, which = parent_and_which bs in
      c_reinit parent which tb0 yb0;
      (match iter_type with
       | Some iter -> set_iter_type bs iter yb0
       | None -> ())


    let get_work_space bs = Cvode.get_work_space (tosession bs)

    let get_num_steps bs = Cvode.get_num_steps (tosession bs)

    let get_num_rhs_evals bs = Cvode.get_num_rhs_evals (tosession bs)

    let get_num_lin_solv_setups bs =
      Cvode.get_num_lin_solv_setups (tosession bs)

    let get_num_err_test_fails bs =
      Cvode.get_num_err_test_fails (tosession bs)

    let get_last_order bs = Cvode.get_last_order (tosession bs)

    let get_current_order bs = Cvode.get_current_order (tosession bs)

    let get_last_step bs = Cvode.get_last_step (tosession bs)

    let get_current_step bs = Cvode.get_current_step (tosession bs)

    let get_actual_init_step bs =
      Cvode.get_actual_init_step (tosession bs)

    let get_current_time bs = Cvode.get_current_time (tosession bs)

    let get_num_stab_lim_order_reds bs =
      Cvode.get_num_stab_lim_order_reds (tosession bs)

    let get_tol_scale_factor bs =
      Cvode.get_tol_scale_factor (tosession bs)

    let get_err_weights bs = Cvode.get_err_weights (tosession bs)
    let get_est_local_errors bs =
      Cvode.get_est_local_errors (tosession bs)

    let get_integrator_stats bs =
      Cvode.get_integrator_stats (tosession bs)

    let print_integrator_stats bs os =
      Cvode.print_integrator_stats (tosession bs) os

    let get_num_nonlin_solv_iters bs =
      Cvode.get_num_nonlin_solv_iters (tosession bs)

    let get_num_nonlin_solv_conv_fails bs =
      Cvode.get_num_nonlin_solv_conv_fails (tosession bs)

    let get_nonlin_solv_stats bs =
      Cvode.get_nonlin_solv_stats (tosession bs)

    module Quadrature =
      struct
        include QuadratureTypes

        external c_quad_initb
            : ('a, 'k) session -> int -> ('a, 'k) nvector -> unit
            = "c_cvodes_adjquad_initb"
        external c_quad_initbs
            : ('a, 'k) session -> int -> ('a, 'k) nvector -> unit
            = "c_cvodes_adjquad_initbs"

        let init bs mf y0 =
          let parent, which = parent_and_which bs in
          let se = bwdsensext bs in
          match mf with
           | NoSens f -> (se.bquadrhsfn <- f;
                             c_quad_initb parent which y0)
           | WithSens f -> (se.bquadrhsfn_sens <- f;
                                c_quad_initbs parent which y0)

        external c_reinit : ('a, 'k) session -> int -> ('a, 'k) nvector -> unit
            = "c_cvodes_adjquad_reinit"

        let reinit bs yqb0 =
          let parent, which = parent_and_which bs in
          c_reinit parent which yqb0

        external c_get : ('a, 'k) session -> int -> ('a, 'k) nvector -> float
            = "c_cvodes_adjquad_get"

        let get bs yqb =
          let parent, which = parent_and_which bs in
          c_get parent which yqb

        type ('a, 'k) tolerance =
            NoStepSizeControl
          | SStolerances of float * float
          | SVtolerances of float * ('a, 'k) nvector

        external set_err_con : ('a, 'k) session -> int -> bool -> unit
            = "c_cvodes_adjquad_set_err_con"

        external sv_tolerances
            : ('a, 'k) session -> int -> float -> ('a, 'k) nvector -> unit
            = "c_cvodes_adjquad_sv_tolerances"

        external ss_tolerances
            : ('a, 'k) session -> int -> float -> float -> unit
            = "c_cvodes_adjquad_ss_tolerances"

        let set_tolerances bs tol =
          let parent, which = parent_and_which bs in
          match tol with
          | NoStepSizeControl -> set_err_con parent which false
          | SStolerances (rel, abs) -> (ss_tolerances parent which rel abs;
                                        set_err_con parent which true)
          | SVtolerances (rel, abs) -> (sv_tolerances parent which rel abs;
                                        set_err_con parent which true)

        let get_num_rhs_evals bs =
          Quadrature.get_num_rhs_evals (tosession bs)

        let get_num_err_test_fails bs =
          Quadrature.get_num_err_test_fails (tosession bs)

        let get_err_weights bs =
          Quadrature.get_err_weights (tosession bs)

        let get_stats bs = Quadrature.get_stats (tosession bs)
      end
  end


(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "c_cvodes_init_module"

let _ =
  c_init_module
    (* Exceptions must be listed in the same order as
       cvodes_exn_index.  *)
    [|Quadrature.QuadNotInitialized;
      Quadrature.QuadRhsFuncFailure;
      Quadrature.FirstQuadRhsFuncFailure;
      Quadrature.RepeatedQuadRhsFuncFailure;
      Quadrature.UnrecoverableQuadRhsFuncFailure;

      Sensitivity.SensNotInitialized;
      Sensitivity.SensRhsFuncFailure;
      Sensitivity.FirstSensRhsFuncFailure;
      Sensitivity.RepeatedSensRhsFuncFailure;
      Sensitivity.UnrecoverableSensRhsFuncFailure;
      Sensitivity.BadSensIdentifier;

      Sensitivity.Quadrature.QuadSensNotInitialized;
      Sensitivity.Quadrature.QuadSensRhsFuncFailure;
      Sensitivity.Quadrature.FirstQuadSensRhsFuncFailure;
      Sensitivity.Quadrature.RepeatedQuadSensRhsFuncFailure;
      Sensitivity.Quadrature.UnrecoverableQuadSensRhsFuncFailure;

      Adjoint.AdjointNotInitialized;
      Adjoint.NoForwardCall;
      Adjoint.ForwardReinitFailure;
      Adjoint.ForwardFailure;
      Adjoint.NoBackwardProblem;
      Adjoint.BadFinalTime;
      Adjoint.BadOutputTime;
    |]
