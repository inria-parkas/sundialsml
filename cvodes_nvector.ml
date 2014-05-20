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

type 'a bsession = {
    parent   : 'a session;
    bsession : 'a session;
    which    : int;
  }

type 'a quadrhsfn = float -> 'a -> 'a -> unit

type 'a sensrhsfn =
    AllAtOnce of (float -> 'a -> 'a -> 'a array -> 'a array -> 'a -> 'a -> unit)
  | OneByOne of (float -> 'a -> 'a -> int -> 'a -> 'a -> 'a -> 'a -> unit)

type 'a quadsensrhsfn =
   float -> 'a nvector -> 'a nvector -> 'a nvector -> 'a nvector array
         -> 'a nvector -> 'a nvector -> unit

type 'a brhsfn =
        BackBasic of (float -> 'a -> 'a -> 'a -> unit)
      | BackWithSens of (float -> 'a -> 'a array -> 'a -> 'a -> unit)

type 'a bquadrhsfn =
        QuadBasic of (float -> 'a -> 'a -> 'a -> unit)
      | QuadWithSens of (float -> 'a -> 'a array -> 'a -> 'a -> unit)

type ('t, 'a) bjacobian_arg =
  {
    jac_t   : float;
    jac_y   : 'a;
    jac_yb  : 'a;
    jac_fyb : 'a;
    jac_tmp : 't
  }

type 'a bprec_solve_arg =
  {
    rvecB   : 'a;
    gammaB : float;
    deltaB : float;
  }

(*
    (* TODO: not for nvectors *)
type 'a bdense_jac_fn =
      ('a triple_tmp, 'a) bjacobian_arg -> Dls.DenseMatrix.t -> unit

 and 'a bband_jac_fn =
      bandrange -> ('a triple_tmp, 'a) bjacobian_arg -> Dls.DenseMatrix.t -> unit

 and bandrange = { mupper : int; mlower : int }
*)

type 'a sensext = {
    (* Quadrature *)
    mutable quadrhsfn     : 'a quadrhsfn;

    (* Forward *)
    mutable senspvals     : Sundials.real_array option;
                            (* keep a reference to prevent garbage collection *)

    mutable sensrhsfn     : (float -> 'a -> 'a -> 'a array
                               -> 'a array -> 'a -> 'a -> unit);
    mutable sensrhsfn1    : (float -> 'a -> 'a -> int -> 'a
                               -> 'a -> 'a -> 'a -> unit);
    mutable quadsensrhsfn : 'a quadsensrhsfn;
  }

let shouldn't_be_called fcn =
  failwith ("internal error in sundials: " ^ fcn ^ " is called")
(*
(* TODO: include in serial vesrsion *)
let dummy_bdense_jac _ _ = shouldn't_be_called "dummy_dense_jac"
let dummy_bband_jac _ _ _ _ = shouldn't_be_called "dummy_band_jac"
*)
let dummy_bprec_setup _ _ _ = shouldn't_be_called "dummy_prec_setup"
let dummy_bprec_solve _ _ _ = shouldn't_be_called "dummy_prec_solve"
let dummy_bjac_times_vec _ _ _ = shouldn't_be_called "dummy_jac_times_vec"

type 'a bsensext = {
    (* Adjoint *)
    mutable brhsfn        : (float -> 'a -> 'a -> 'a -> unit);
    mutable brhsfn1       : (float -> 'a -> 'a array -> 'a -> 'a -> unit);
    mutable bquadrhsfn    : (float -> 'a -> 'a -> 'a -> unit);
    mutable bquadrhsfn1   : (float -> 'a -> 'a array -> 'a -> 'a -> unit);

    mutable bpresetupfn : ('a triple_tmp, 'a) bjacobian_arg -> bool
                            -> float -> bool;
    mutable bpresolvefn : ('a single_tmp, 'a) bjacobian_arg
                            -> 'a bprec_solve_arg -> 'a -> unit;
    mutable bjactimesfn : ('a single_tmp, 'a) bjacobian_arg -> 'a -> 'a -> unit;

    (*
    (* TODO: include in serial version *)
    mutable bjacfn      : 'a bdense_jac_fn;
    mutable bbandjacfn  : 'a bband_jac_fn;
    *)
  }

exception NoSensExt

let sensext (s : 'a session) =
  match s.sensext with
  | None -> raise NoSensExt
  | Some se -> ((Obj.obj se) : 'a sensext)

let new_sensext s =
  try sensext s
  with NoSensExt ->
    let se = {
        senspvals     = None;
        quadrhsfn     = (fun _ _ _ -> ());
        sensrhsfn     = (fun _ _ _ _ _ _ _ -> ());
        sensrhsfn1    = (fun _ _ _ _ _ _ _ _ -> ());
        quadsensrhsfn = (fun _ _ _ _ _ _ _ -> ());
      }
    in
    s.sensext <- Some (Obj.repr se);
    se

let bsensext {bsession=bs} =
  match bs.sensext with
  | None -> raise NoSensExt
  | Some se -> ((Obj.obj se) : 'a bsensext)

let empty_bsensext = {
      brhsfn      = (fun _ _ _ _ -> ());
      brhsfn1     = (fun _ _ _ _ _ -> ());
      bquadrhsfn  = (fun _ _ _ _ -> ());
      bquadrhsfn1 = (fun _ _ _ _ _ -> ());

      bpresetupfn = dummy_bprec_setup;
      bpresolvefn = dummy_bprec_solve;
      bjactimesfn = dummy_bjac_times_vec;

(* TODO: for serial
      bdense_jac = dummy_bdense_jac;
      bband_jac = dummy_bband_jac;
*)
    }

let new_bsensext ({bsession=bs} as s) =
  try bsensext s
  with NoSensExt ->
    let se = empty_bsensext in
    bs.sensext <- Some (Obj.repr se);
    se

(* TODO: add callback 'trampolines' *)

module Quadrature =
  struct

    (* TODO: add to the standard CHECK_FLAG function? *)
    exception QuadNotInitialized
    exception QuadRhsFuncFailure
    exception FirstQuadRhsFuncErr
    exception RepeatedQuadRhsFuncErr
    exception UnrecoverableQuadRhsFuncErr

    type 'a quadrhsfn = float -> 'a -> 'a -> unit

    external c_init : 'a session -> 'a nvector -> unit
        = "c_nvec_cvodes_quad_init"

    let init s f v0 =
      let se = new_sensext s in
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

module Forward =
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
    

    (* TODO: add to the standard CHECK_FLAG function? *)
    
    exception SensNotInitialized
    exception SensRhsFuncFailure
    exception FirstSensRhsFuncErr
    exception RepeatedSensRhsFuncErr
    exception UnrecoverableSensRhsFuncErr

    type fwd_method =
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

    external c_init : 'a session -> fwd_method -> 'a nvector array -> unit
        = "c_nvec_cvodes_sens_init"

    external c_init_1 : 'a session -> fwd_method -> 'a nvector array -> unit
        = "c_nvec_cvodes_sens_init_1"

    (* TODO: check that pbar and plist are ns long. *)
    external set_params : 'a session -> sens_params -> unit
        = "c_cvodes_sens_set_params"

    let init s tol fmethod sparams fm v0 =
      let se = sensext s in
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

    external reinit : 'a session -> fwd_method -> 'a nvector array -> unit
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

        (* TODO: add to the standard CHECK_FLAG function? *)
        exception QuadSensNotInitialized
        exception QuadSensRhsFuncFailure
        exception FirstQuadSensRhsFuncErr
        exception RepeatedQuadSensRhsFuncErr
        exception UnrecoverableQuadSensRhsFuncErr

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
          let se = sensext s in
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

    (* TODO: add to the standard CHECK_FLAG function? *)
    (* TODO: check that all return flags have been considered; here and
       elsewhere *)

    exception SensNotInitialized
    exception NoForwardCall
    exception ForwardReinitializationFailed
    exception ForwardFailed
    exception NoBackwardProblem
    exception NoBackwardProblem

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

    type 'a tolerance =
      | SSTolerances of float * float
      | SVTolerances of float * 'a nvector

    external ss_tolerances
        : 'a session -> int -> float -> float -> unit
        = "c_cvodes_adj_ss_tolerancesb"

    external sv_tolerances
        : 'a session -> int -> float -> 'a nvector -> unit
        = "c_nvec_cvodes_adj_sv_tolerancesb"

    let set_tolerances {parent=s; which=w} tol =
      match tol with
      | SSTolerances (rel, abs) -> ss_tolerances s w rel abs
      | SVTolerances (rel, abs) -> sv_tolerances s w rel abs

    type 'a _bsession = 'a bsession
    type 'a bsession = 'a _bsession

    external c_diag : 'a session -> int -> unit
      = "c_cvode_adj_diag"

    external c_spils_set_preconditioner
      : 'a session -> int -> bool -> bool -> unit
      = "c_nvec_cvode_adj_spils_set_preconditioner"

    external c_spils_spgmr
      : 'a session -> int -> int -> Spils.preconditioning_type -> unit
      = "c_cvode_adj_spils_spgmr"

    external c_spils_spbcg
      : 'a session -> int -> int -> Spils.preconditioning_type -> unit
      = "c_cvode_adj_spils_spbcg"

    external c_spils_sptfqmr
      : 'a session -> int -> int -> Spils.preconditioning_type -> unit
      = "c_cvode_adj_spils_sptfqmr"

(* XXX serial only
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

    external c_set_functional : 'a session -> int -> unit
      = "c_cvode_adj_set_functional"

    (* TODO: rework for serial nvectors *)
    let set_iter_type ({parent=s;which=w} as session) iter =
      let se = bsensext session in
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
            c_spils_set_preconditioner s w
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
      | Functional -> c_set_functional s w
      | Newton linsolv ->
        (* Iter type will be set to CV_NEWTON in the functions that set the linear
           solver.  *)
        match linsolv with
        | Diag -> c_diag s w
        | Spgmr (par, cb) ->
            let maxl = match par.maxl with None -> 0 | Some ml -> ml in
            c_spils_spgmr s w maxl par.prec_type;
            set_precond par.prec_type cb
        | Spbcg (par, cb) ->
            let maxl = match par.maxl with None -> 0 | Some ml -> ml in
            c_spils_spbcg s w maxl par.prec_type;
            set_precond par.prec_type cb
        | Sptfqmr (par, cb) ->
            let maxl = match par.maxl with None -> 0 | Some ml -> ml in
            c_spils_sptfqmr s w maxl par.prec_type;
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
      let se, (cvode_mem, which, backref, err_file) =
        match mf with
        | BackBasic f -> begin
              { empty_bsensext with brhsfn = f },
              c_init_backward weakref lmm iter y0 t0
            end
        | BackWithSens f -> begin
              { empty_bsensext with brhsfn1 = f },
              c_init_backward_1 weakref lmm iter y0 t0
            end
      in
      (* cvode_mem and backref have to be immediately captured in a session and
         associated with the finalizer before we do anything else.  *)
      let backsession = {
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

              sensext    = Some (Obj.repr se);
            } in
      Gc.finalise bsession_finalize backsession;
      Weak.set weakref 0 (Some backsession);
      (* Now the session is safe to use.  If any of the following fails and raises
         an exception, the GC will take care of freeing cvode_mem and backref.  *)
      let session = {parent=s; bsession=backsession; which=which} in
      set_iter_type session iter;
      set_tolerances session tol;
      session

    external c_reinit : 'a session -> int -> float -> 'a nvector -> unit
        = "c_nvec_cvodes_adj_reinit"
    let reinit {parent=s; which=w} = c_reinit s w

    external c_backward_normal : 'a session -> float -> unit
        = "c_cvodes_adj_backward_normal"
    let backward_normal {parent=s} = c_backward_normal s

    external c_backward_one_step : 'a session -> float -> unit
        = "c_cvodes_adj_backward_one_step"
    let backward_one_step {parent=s} = c_backward_one_step s

    external c_get : 'a session -> int -> 'a nvector -> float
        = "c_nvec_cvodes_adj_backward_one_step"
    let get {parent=s; which=w} = c_get s w

    let get_dky {bsession=bs} = Cvode_nvector.get_dky bs

    external set_no_sensitivity : 'a session -> unit
        = "c_cvodes_adj_set_no_sensitivity"

    external c_set_max_ord : 'a session -> int -> int -> unit
        = "c_cvodes_set_max_ord"
    let set_max_ord {parent=s; which=w} = c_set_max_ord s w

    external c_set_max_num_steps : 'a session -> int -> int -> unit
        = "c_cvodes_set_max_num_steps"
    let set_max_num_steps {parent=s; which=w} = c_set_max_num_steps s w

    external c_set_init_step : 'a session -> int -> float -> unit
        = "c_cvodes_set_init_step"
    let set_init_step {parent=s; which=w} = c_set_init_step s w

    external c_set_min_step : 'a session -> int -> float -> unit
        = "c_cvodes_set_min_step"
    let set_min_step {parent=s; which=w} = c_set_min_step s w

    external c_set_max_step : 'a session -> int -> float -> unit
        = "c_cvodes_set_max_step"
    let set_max_step {parent=s; which=w} = c_set_max_step s w

    external c_set_stab_lim_det : 'a session -> int -> bool -> unit
        = "c_cvodes_set_stab_lim_det"
    let set_stab_lim_det {parent=s; which=w} = c_set_stab_lim_det s w

    module Diag =
      struct
        let get_work_space {bsession=s} =
          Cvode_nvector.Diag.get_work_space s

        let get_num_rhs_evals {bsession=s} =
          Cvode_nvector.Diag.get_num_rhs_evals s
      end

    module Spils =
      struct
        external c_set_prec_type
            : 'a session -> int -> Spils.preconditioning_type -> unit
            = "c_cvodes_set_prec_type"
        let set_prec_type {parent=s; which=w} = c_set_prec_type s w

        external c_set_gs_type
            : 'a session -> int -> Spils.gramschmidt_type -> unit
            = "c_cvodes_set_gs_type"
        let set_gs_type {parent=s; which=w} = c_set_gs_type s w

        external c_set_eps_lin : 'a session -> int -> float -> unit
            = "c_cvodes_set_eps_lin"
        let set_eps_lin {parent=s; which=w} = c_set_eps_lin s w

        external c_set_maxl : 'a session -> int -> int -> unit
            = "c_cvodes_set_maxl"
        let set_maxl {parent=s; which=w} omaxl =
          let maxl = match omaxl with
                     | None -> 0
                     | Some x -> x
          in c_set_maxl s w maxl

        let get_work_space {bsession=s} =
          Cvode_nvector.Spils.get_work_space s

        let get_num_lin_iters {bsession=s} =
          Cvode_nvector.Spils.get_num_lin_iters s

        let get_num_conv_fails {bsession=s} =
          Cvode_nvector.Spils.get_num_conv_fails s

        let get_num_prec_evals {bsession=s} =
          Cvode_nvector.Spils.get_num_prec_evals s

        let get_num_prec_solves {bsession=s} =
          Cvode_nvector.Spils.get_num_prec_solves s

        let get_num_jtimes_evals {bsession=s} =
          Cvode_nvector.Spils.get_num_jtimes_evals s

        let get_num_rhs_evals {bsession=s} =
          Cvode_nvector.Spils.get_num_rhs_evals s
      end

    module BandPrec =
      struct
        let get_work_space {bsession=s} =
          Cvode_nvector.BandPrec.get_work_space s
        let get_num_rhs_evals {bsession=s} =
          Cvode_nvector.BandPrec.get_num_rhs_evals s
      end

    let get_work_space {bsession=s} = Cvode_nvector.get_work_space s

    let get_num_steps {bsession=s} = Cvode_nvector.get_num_steps s

    let get_num_rhs_evals {bsession=s} = Cvode_nvector.get_num_rhs_evals s

    let get_num_lin_solv_setups {bsession=s} =
      Cvode_nvector.get_num_lin_solv_setups s

    let get_num_err_test_fails {bsession=s} =
      Cvode_nvector.get_num_err_test_fails s

    let get_last_order {bsession=s} = Cvode_nvector.get_last_order s

    let get_current_order {bsession=s} = Cvode_nvector.get_current_order s

    let get_last_step {bsession=s} = Cvode_nvector.get_last_step s

    let get_current_step {bsession=s} = Cvode_nvector.get_current_step s

    let get_actual_init_step {bsession=s} = Cvode_nvector.get_actual_init_step s

    let get_current_time {bsession=s} = Cvode_nvector.get_current_time s

    let get_num_stab_lim_order_reds {bsession=s} =
      Cvode_nvector.get_num_stab_lim_order_reds s

    let get_tol_scale_factor {bsession=s} = Cvode_nvector.get_tol_scale_factor s

    let get_err_weights {bsession=s} = Cvode_nvector.get_err_weights s
    let get_est_local_errors {bsession=s} = Cvode_nvector.get_est_local_errors s

    let get_integrator_stats {bsession=s} = Cvode_nvector.get_integrator_stats s

    let print_integrator_stats {bsession=s} =
      Cvode_nvector.print_integrator_stats s

    let get_num_nonlin_solv_iters {bsession=s} =
      Cvode_nvector.get_num_nonlin_solv_iters s

    let get_num_nonlin_solv_conv_fails {bsession=s} =
      Cvode_nvector.get_num_nonlin_solv_conv_fails s

    module Quadrature =
      struct
        type 'a _bquadrhsfn = 'a bquadrhsfn =
            QuadBasic of (float -> 'a -> 'a -> 'a -> unit)
          | QuadWithSens of (float -> 'a -> 'a array -> 'a -> 'a -> unit)
        type 'a bquadrhsfn = 'a _bquadrhsfn =
            QuadBasic of (float -> 'a -> 'a -> 'a -> unit)
          | QuadWithSens of (float -> 'a -> 'a array -> 'a -> 'a -> unit)

        external c_init : 'a session -> int -> 'a nvector -> unit
            = "c_nvec_cvodes_quad_initb"
        external c_init_s : 'a session -> int -> 'a nvector -> unit
            = "c_nvec_cvodes_quad_initbs"

        let init ({parent=p; which=w} as s) mf y0 =
          let bse = new_bsensext s in
          match mf with
           | QuadBasic f -> (bse.bquadrhsfn <- f;
                             c_init p w y0)
           | QuadWithSens f -> (bse.bquadrhsfn1 <- f;
                                c_init_s p w y0)

        external c_reinit : 'a session -> int -> 'a nvector -> unit
            = "c_nvec_cvodes_quad_reinitb"

        let reinit {parent=s; which=w} = c_reinit s w

        external c_get : 'a session -> int -> 'a nvector array -> float
            = "c_nvec_cvodes_quad_getb"

        let get {parent=s; which=w} = c_get s w

        type 'a tolerance =
            NoStepSizeControl
          | SSTolerances of float * float
          | SVTolerances of float * 'a nvector

        external set_err_con : 'a session -> int -> bool -> unit
            = "c_nvec_cvodes_quad_set_err_conb"

        external sv_tolerances
            : 'a session -> int -> float -> 'a nvector -> unit
            = "c_nvec_cvodes_quad_sv_tolerancesb"

        external ss_tolerances  : 'a session -> int -> float -> float -> unit
            = "c_cvodes_quad_ss_tolerancesb"

        let set_tolerances {parent=s; which=w} tol =
          match tol with
          | NoStepSizeControl -> set_err_con s w false
          | SSTolerances (rel, abs) -> (set_err_con s w true;
                                        ss_tolerances s w rel abs)
          | SVTolerances (rel, abs) -> (set_err_con s w true;
                                        sv_tolerances s w rel abs)

        let get_num_rhs_evals {bsession=s} =
          Quadrature.get_num_rhs_evals s
        let get_num_err_test_fails {bsession=s} =
          Quadrature.get_num_err_test_fails s
        let get_quad_err_weights {bsession=s} =
          Quadrature.get_quad_err_weights s
        let get_stats {bsession=s} = Quadrature.get_stats s
      end
  end

