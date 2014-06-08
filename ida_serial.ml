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

include Ida

type nvec = Sundials.Carray.t
type val_array = Sundials.Carray.t
type der_array = Sundials.Carray.t

type root_array = Sundials.Roots.t
type root_val_array = Sundials.Roots.val_array

type single_tmp = val_array
type double_tmp = val_array * val_array
type triple_tmp = val_array * val_array * val_array

type 't jacobian_arg =
  {
    jac_t    : float;
    jac_y    : val_array;
    jac_y'   : der_array;
    jac_res  : val_array;
    jac_coef : float;
    jac_tmp  : 't
  }

(* Note this definition differs from the one in ida_nvector, so the
   tag values are different.  *)
type linear_solver =
  | Dense of dense_jac_fn option
  | LapackDense of dense_jac_fn option
  | Band of bandrange * band_jac_fn option
  | LapackBand of bandrange * band_jac_fn option
  | Spgmr of spils_params
  | Spbcg of spils_params
  | Sptfqmr of spils_params
and dense_jac_fn = triple_tmp jacobian_arg -> Dls.DenseMatrix.t -> unit
and band_jac_fn = bandrange -> triple_tmp jacobian_arg
                            -> Dls.BandMatrix.t -> unit
and spils_params =
  {
    maxl : int option;
    prec_solve_fn : (single_tmp jacobian_arg -> val_array -> val_array -> float
                     -> unit) option;
    prec_setup_fn : (triple_tmp jacobian_arg -> unit) option;
    jac_times_vec_fn : (double_tmp jacobian_arg -> val_array -> val_array
                        -> unit) option;
  }

let spils_no_precond = { maxl = None;
                         prec_solve_fn = None;
                         prec_setup_fn = None;
                         jac_times_vec_fn = None;
                       }

type ida_mem
type c_weak_ref
type ida_file
type session = {
        ida        : ida_mem;
        backref    : c_weak_ref; (* Back reference to the whole record that
                                    goes through a generational root and a
                                    weak pointer, for access from callbacks *)
        neqs       : int;
        nroots     : int;
        err_file   : ida_file;

        (* Temporary storage for exceptions raised within callbacks.  *)
        mutable exn_temp   : exn option;

        mutable resfn      : float -> val_array -> der_array -> val_array
                             -> unit;
        mutable rootsfn    : float -> val_array -> der_array -> root_val_array
                             -> unit;
        mutable errh       : Sundials.error_details -> unit;
        mutable errw       : val_array -> nvec -> unit;
        mutable jacfn      : dense_jac_fn;
        mutable bandjacfn  : band_jac_fn;
        mutable presetupfn : triple_tmp jacobian_arg -> unit;
        mutable presolvefn : single_tmp jacobian_arg -> val_array -> val_array
                               -> float -> unit;
        mutable jactimesfn : double_tmp jacobian_arg -> val_array -> val_array
                               -> unit;

        (* To be manipulated from the C side only.  *)
        mutable safety_check_flags : int;
      }

external sv_tolerances  : session -> float -> nvec -> unit
  = "c_ba_ida_sv_tolerances"
external ss_tolerances  : session -> float -> float -> unit
  = "c_ida_ss_tolerances"
external wf_tolerances  : session -> unit
  = "c_ba_ida_wf_tolerances"

type tolerance =
  | SSTolerances of float * float
    (** [(rel, abs)] : scalar relative and absolute tolerances. *)
  | SVTolerances of float * nvec
    (** [(rel, abs)] : scalar relative and vector absolute tolerances. *)
  | WFTolerances of (val_array -> val_array -> unit)
    (** Specifies a function [efun y ewt] that sets the multiplicative
        error weights Wi for use in the weighted RMS norm. The function is
        passed the dependent variable vector [y] and is expected to set the
        values inside the error-weight vector [ewt]. *)

let default_tolerances = SSTolerances (1.0e-4, 1.0e-8)

let set_tolerances s tol =
  match tol with
  | SSTolerances (rel, abs) -> ss_tolerances s rel abs
  | SVTolerances (rel, abs) -> sv_tolerances s rel abs
  | WFTolerances ferrw -> (s.errw <- ferrw; wf_tolerances s)

let read_weak_ref x : session =
  match Weak.get x 0 with
  | Some y -> y
  | None -> raise (Failure "Internal error: weak reference is dead")
let adjust_retcode = fun session check_recoverable f x ->
  try f x; 0
  with
  | Sundials.RecoverableFailure when check_recoverable -> 1
  | e -> (session.exn_temp <- Some e; -1)
let call_resfn session t y y' res =
  let session = read_weak_ref session in
  adjust_retcode session true (session.resfn t y y') res
let call_rootsfn session t y y' roots =
  let session = read_weak_ref session in
  adjust_retcode session false (session.rootsfn t y y') roots
let call_errw session y ewt =
  let session = read_weak_ref session in
  adjust_retcode session false (session.errw y) ewt
let call_errh session details =
  let session = read_weak_ref session in
  try session.errh details
  with e ->
    prerr_endline ("Warning: error handler function raised an exception.  " ^
                   "This exception will not be propagated.")
let call_jacfn session jac j =
  let session = read_weak_ref session in
  adjust_retcode session true (session.jacfn jac) j
let call_bandjacfn session range jac j =
  let session = read_weak_ref session in
  adjust_retcode session true (session.bandjacfn range jac) j
let call_presetupfn session jac =
  let session = read_weak_ref session in
  adjust_retcode session true session.presetupfn jac
let call_presolvefn session jac r z delta =
  let session = read_weak_ref session in
  adjust_retcode session true (session.presolvefn jac r z) delta
let call_jactimesfn session jac r z =
  let session = read_weak_ref session in
  adjust_retcode session true (session.jactimesfn jac r) z
let _ =
  Callback.register "c_ba_ida_call_resfn"         call_resfn;
  Callback.register "c_ba_ida_call_rootsfn"       call_rootsfn;
  Callback.register "c_ba_ida_call_errh"          call_errh;
  Callback.register "c_ba_ida_call_errw"          call_errw;
  Callback.register "c_ba_ida_call_jacfn"         call_jacfn;
  Callback.register "c_ba_ida_call_bandjacfn"     call_bandjacfn;
  Callback.register "c_ba_ida_call_presetupfn"    call_presetupfn;
  Callback.register "c_ba_ida_call_presolvefn"    call_presolvefn;
  Callback.register "c_ba_ida_call_jactimesfn"    call_jactimesfn;

external c_init
  : 'a Weak.t -> float -> val_array -> der_array -> (ida_mem * c_weak_ref * ida_file)
  = "c_ba_ida_init"
external session_finalize : session -> unit
  = "c_ida_session_finalize"

external c_root_init : session -> int -> unit
  = "c_ba_ida_root_init"
let root_init ida (nroots, rootsfn) =
  c_root_init ida nroots;
  ida.rootsfn <- rootsfn

external c_dls_dense : session -> bool -> unit
  = "c_ba_ida_dls_dense"

external c_dls_lapack_dense : session -> bool -> unit
  = "c_ba_ida_dls_lapack_dense"

external c_dls_band : session -> int -> int -> bool -> unit
  = "c_ba_ida_dls_band"

external c_dls_lapack_band : session -> int -> int -> bool -> unit
  = "c_ba_ida_dls_lapack_band"

external c_spils_spgmr : session -> int -> unit
  = "c_ba_ida_spils_spgmr"

external c_spils_spbcg : session -> int -> unit
  = "c_ba_ida_spils_spbcg"

external c_spils_sptfqmr : session -> int -> unit
  = "c_ba_ida_spils_sptfqmr"

external c_spils_set_preconditioner : session -> bool -> bool -> unit
  = "c_ba_ida_spils_set_preconditioner"

let shouldn't_be_called fcn =
  failwith ("internal error in sundials: " ^ fcn ^ " is called")
let dummy_dense_jac _ _ = shouldn't_be_called "dummy_dense_jac"
let dummy_band_jac _ _ _ = shouldn't_be_called "dummy_band_jac"
let dummy_prec_setup _ = shouldn't_be_called "dummy_prec_setup"
let dummy_prec_solve _ _ _ _ = shouldn't_be_called "dummy_prec_solve"
let dummy_jac_times_vec _ _ _ = shouldn't_be_called "dummy_jac_times_vec"

let set_linear_solver session solver =
  let optionally f = function
    | None -> ()
    | Some x -> f x
  in
  let set_precond p =
    match p.prec_solve_fn with
    | None -> ()
    | Some solve_fn ->
      c_spils_set_preconditioner session
        (p.prec_setup_fn <> None)
        (p.jac_times_vec_fn <> None);
      session.presolvefn <- solve_fn;
      optionally (fun f -> session.presetupfn <- f) p.prec_setup_fn;
      optionally (fun f -> session.jactimesfn <- f) p.jac_times_vec_fn
  in
  (* Release references to all linear-solver--related callbacks.  *)
  session.jacfn      <- dummy_dense_jac;
  session.bandjacfn  <- dummy_band_jac;
  session.presetupfn <- dummy_prec_setup;
  session.presolvefn <- dummy_prec_solve;
  session.jactimesfn <- dummy_jac_times_vec;
  match solver with
  | Dense jac ->
    c_dls_dense session (jac <> None);
    optionally (fun f -> session.jacfn <- f) jac
  | LapackDense jac ->
    c_dls_lapack_dense session (jac <> None);
    optionally (fun f -> session.jacfn <- f) jac
  | Band (p, jac) ->
    c_dls_band session p.mupper p.mlower (jac <> None);
    optionally (fun f -> session.bandjacfn <- f) jac
  | LapackBand (p, jac) ->
    c_dls_lapack_band session p.mupper p.mlower (jac <> None);
    optionally (fun f -> session.bandjacfn <- f) jac
  | Spgmr p ->
      let maxl = match p.maxl with None -> 0 | Some ml -> ml in
      c_spils_spgmr session maxl;
      set_precond p
  | Spbcg p ->
      let maxl = match p.maxl with None -> 0 | Some ml -> ml in
      c_spils_spbcg session maxl;
      set_precond p
  | Sptfqmr p ->
      let maxl = match p.maxl with None -> 0 | Some ml -> ml in
      c_spils_sptfqmr session maxl;
      set_precond p

let init linsolv tol resfn ?(roots=no_roots) ?(t0=0.) y y' =
  let (nroots, rootsfn) = roots in
  if nroots < 0 then
    raise (Invalid_argument "number of root functions is negative");
  let neqs = Sundials.Carray.length y in
  (* IDA doesn't check if y and y' have the same length, and corrupt memory if
   * they don't.  *)
  if neqs <> Sundials.Carray.length y' then
    raise (Invalid_argument "y and y' have inconsistent sizes");
  let weakref = Weak.create 1 in
  let (ida_mem, backref, err_file) = c_init weakref t0 y y' in
  (* ida_mem and backref have to be immediately captured in a session and
     associated with the finalizer before we do anything else.  *)
  let session = { ida        = ida_mem;
                  backref    = backref;
                  neqs       = neqs;
                  nroots     = nroots;
                  err_file   = err_file;
                  exn_temp   = None;
                  resfn      = resfn;
                  rootsfn    = rootsfn;
                  errh       = (fun _ -> ());
                  errw       = (fun _ _ -> ());
                  jacfn      = dummy_dense_jac;
                  bandjacfn  = dummy_band_jac;
                  presetupfn = dummy_prec_setup;
                  presolvefn = dummy_prec_solve;
                  jactimesfn = dummy_jac_times_vec;
                  safety_check_flags = 0;
                }
  in
  Gc.finalise session_finalize session;
  Weak.set weakref 0 (Some session);
  (* Now the session is safe to use.  If any of the following fails and raises
     an exception, the GC will take care of freeing ida_mem and backref.  *)
  if nroots > 0 then
    c_root_init session nroots;
  set_linear_solver session linsolv;
  set_tolerances session tol;
  session

let nroots { nroots } = nroots
let neqs { neqs } = neqs

external c_reinit : session -> float -> val_array -> der_array -> unit
    = "c_ba_ida_reinit"

let reinit ida ?linsolv ?roots t0 y0 y'0 =
  c_reinit ida t0 y0 y'0;
  (match linsolv with
   | None -> ()
   | Some linsolv -> set_linear_solver ida linsolv);
  (match roots with
   | None -> ()
   | Some roots -> root_init ida roots)

external get_root_info  : session -> root_array -> unit
    = "c_ida_get_root_info"

external solve_normal
    : session -> float -> val_array -> der_array -> float * solver_result
    = "c_ba_ida_normal"

external solve_one_step
    : session -> float -> val_array -> der_array -> float * solver_result
    = "c_ba_ida_one_step"

external get_dky
    : session -> float -> int -> nvec -> unit
    = "c_ba_ida_get_dky"

external get_integrator_stats   : session -> integrator_stats
    = "c_ida_get_integrator_stats"

external get_work_space         : session -> int * int
    = "c_ida_get_work_space"

external get_num_steps          : session -> int
    = "c_ida_get_num_steps"

external get_num_res_evals      : session -> int
    = "c_ida_get_num_res_evals"

external get_num_lin_solv_setups : session -> int
    = "c_ida_get_num_lin_solv_setups"

external get_num_err_test_fails : session -> int
    = "c_ida_get_num_err_test_fails"

external get_last_order         : session -> int
    = "c_ida_get_last_order"

external get_current_order      : session -> int
    = "c_ida_get_current_order"

external get_actual_init_step   : session -> float
    = "c_ida_get_actual_init_step"

external get_last_step          : session -> float
    = "c_ida_get_last_step"

external get_current_step       : session -> float
    = "c_ida_get_current_step"

external get_current_time       : session -> float
    = "c_ida_get_current_time"

let print_integrator_stats s =
  let stats = get_integrator_stats s
  in
    Printf.printf "num_steps = %d\n"           stats.num_steps;
    Printf.printf "num_res_evals = %d\n"       stats.num_res_evals;
    Printf.printf "num_lin_solv_setups = %d\n" stats.num_lin_solv_setups;
    Printf.printf "num_err_test_fails = %d\n"  stats.num_err_test_fails;
    Printf.printf "last_order = %d\n"          stats.last_order;
    Printf.printf "current_order = %d\n"       stats.current_order;
    Printf.printf "actual_init_step = %e\n"    stats.actual_init_step;
    Printf.printf "last_step = %e\n"           stats.last_step;
    Printf.printf "current_step = %e\n"        stats.current_step;
    Printf.printf "current_time = %e\n"        stats.current_time;

external set_error_file : session -> string -> bool -> unit
    = "c_ida_set_error_file"

external set_err_handler_fn  : session -> unit
    = "c_ba_ida_set_err_handler_fn"

let set_err_handler_fn s ferrh =
  s.errh <- ferrh;
  set_err_handler_fn s

external clear_err_handler_fn  : session -> unit
    = "c_ba_ida_clear_err_handler_fn"

let clear_err_handler_fn s =
  s.errh <- (fun _ -> ());
  clear_err_handler_fn s

external set_max_ord            : session -> int -> unit
    = "c_ida_set_max_ord"
external set_max_num_steps      : session -> int -> unit
    = "c_ida_set_max_num_steps"
external set_init_step          : session -> float -> unit
    = "c_ida_set_init_step"
external set_max_step           : session -> float -> unit
    = "c_ida_set_max_step"
external set_stop_time          : session -> float -> unit
    = "c_ida_set_stop_time"
external set_max_err_test_fails : session -> int -> unit
    = "c_ida_set_max_err_test_fails"
external set_max_nonlin_iters   : session -> int -> unit
    = "c_ida_set_max_nonlin_iters"
external set_max_conv_fails     : session -> int -> unit
    = "c_ida_set_max_conv_fails"
external set_nonlin_conv_coef   : session -> float -> unit
    = "c_ida_set_nonlin_conv_coef"

external set_root_direction'    : session -> RootDirs.t -> unit
    = "c_ida_set_root_direction"

let set_root_direction s rda = 
  set_root_direction' s (RootDirs.copy_n (nroots s) rda)

let set_all_root_directions s rd =
  set_root_direction' s (RootDirs.make (nroots s) rd)

external set_no_inactive_root_warn      : session -> unit
    = "c_ida_set_no_inactive_root_warn"

external get_tol_scale_factor           : session -> float
    = "c_ida_get_tol_scale_factor"

external get_err_weights                : session -> nvec -> unit
    = "c_ba_ida_get_err_weights"

external get_est_local_errors           : session -> nvec -> unit
    = "c_ba_ida_get_est_local_errors"

external get_num_nonlin_solv_iters      : session -> int
    = "c_ida_get_num_nonlin_solv_iters"

external get_num_nonlin_solv_conv_fails : session -> int
    = "c_ida_get_num_nonlin_solv_conv_fails"

external get_num_g_evals                : session -> int
    = "c_ida_get_num_g_evals"

module Dls =
  struct

    external set_dense_jac_fn  : session -> unit
        = "c_ba_ida_dls_set_dense_jac_fn"

    let set_dense_jac_fn s fjacfn =
      s.jacfn <- fjacfn;
      set_dense_jac_fn s

    external clear_dense_jac_fn : session -> unit
        = "c_ba_ida_dls_clear_dense_jac_fn"

    let clear_dense_jac_fn s =
      s.jacfn <- dummy_dense_jac;
      clear_dense_jac_fn s

    external set_band_jac_fn   : session -> unit
        = "c_ba_ida_dls_set_band_jac_fn"

    let set_band_jac_fn s fbandjacfn =
      s.bandjacfn <- fbandjacfn;
      set_band_jac_fn s

    external clear_band_jac_fn : session -> unit
        = "c_ba_ida_dls_clear_band_jac_fn"

    let clear_band_jac_fn s =
      s.bandjacfn <- dummy_band_jac;
      clear_band_jac_fn s

    external get_work_space : session -> int * int
        = "c_ida_dls_get_work_space"

    external get_num_jac_evals    : session -> int
        = "c_ida_dls_get_num_jac_evals"

    external get_num_res_evals    : session -> int
        = "c_ida_dls_get_num_res_evals"
  end

module Spils =
  struct
    external set_preconditioner  : session -> unit
        = "c_ba_ida_set_preconditioner"

    let set_preconditioner s fpresetupfn fpresolvefn =
      s.presetupfn <- fpresetupfn;
      s.presolvefn <- fpresolvefn;
      set_preconditioner s

    external set_jac_times_vec_fn : session -> unit
        = "c_ba_ida_set_jac_times_vec_fn"

    let set_jac_times_vec_fn s fjactimesfn =
      s.jactimesfn <- fjactimesfn;
      set_jac_times_vec_fn s

    external clear_jac_times_vec_fn : session -> unit
        = "c_ba_ida_clear_jac_times_vec_fn"

    let clear_jac_times_vec_fn s =
      s.jactimesfn <- (fun _ _ _ -> ());
      clear_jac_times_vec_fn s

    external set_gs_type : session -> Spils.gramschmidt_type -> unit
        = "c_ida_set_gs_type"

    external set_eps_lin            : session -> float -> unit
        = "c_ida_set_eps_lin"

    external set_maxl               : session -> int -> unit
        = "c_ida_set_maxl"

    external get_num_lin_iters      : session -> int
        = "c_ida_spils_get_num_lin_iters"

    external get_num_conv_fails     : session -> int
        = "c_ida_spils_get_num_conv_fails"

    external get_work_space         : session -> int * int
        = "c_ida_spils_get_work_space"

    external get_num_prec_evals     : session -> int
        = "c_ida_spils_get_num_prec_evals"

    external get_num_prec_solves    : session -> int
        = "c_ida_spils_get_num_prec_solves"

    external get_num_jtimes_evals   : session -> int
        = "c_ida_spils_get_num_jtimes_evals"

    external get_num_res_evals      : session -> int
        = "c_ida_spils_get_num_res_evals"

  end

module Constraints =
  struct
    type t = nvec
    type constraint_type =
    | Unconstrained
    | NonNegative
    | NonPositive
    | Positive
    | Negative

    let constraint_type_of_float = function
      | 0.0  -> Unconstrained
      | 1.0  -> NonNegative
      | -1.0 -> NonPositive
      | 2.0  -> Positive
      | -2.0 -> Negative
      | f -> raise (Invalid_argument
                      ("invalid constraint: " ^ string_of_float f))
    let float_of_constraint_type = function
      | Unconstrained -> 0.0
      | NonNegative   -> 1.0
      | NonPositive   -> -1.0
      | Positive      -> 2.0
      | Negative      -> -2.0

    let create = Carray.create
    let init n v = Carray.init n (float_of_constraint_type v)
    let length = Carray.length
    let of_array a =
      let ret = create (Array.length a) in
      for i = 0 to Array.length a - 1 do
        ret.{i} <- float_of_constraint_type a.(i)
      done;
      ret

    let get a i = constraint_type_of_float a.{i}
    let set a i x = a.{i} <- float_of_constraint_type x
    let fill a t =
      let x = float_of_constraint_type t in
      Carray.fill a x

    let blit a b = Carray.blit a b
  end

external set_constraints : session -> Constraints.t -> unit
  = "c_ba_ida_set_constraints"

external set_id : session -> Id.t -> unit
  = "c_ba_ida_set_id"
let set_var_types = set_id
external set_suppress_alg : session -> bool -> unit
  = "c_ida_set_suppress_alg"

external get_num_backtrack_ops : session -> int
  = "c_ida_get_num_backtrack_ops"

external c_calc_ic_y : session -> val_array option -> float -> unit
  = "c_ba_ida_calc_ic_y"
let calc_ic_y session ?y tout1 =
  let len = function
    | None -> neqs session
    | Some x -> Carray.length x
  in
  if len y <> neqs session then
    invalid_arg "calc_ic_y: ~y has incorrect length";
  c_calc_ic_y session y tout1

external c_calc_ic_ya_yd' :
  session -> val_array option -> der_array option -> Id.t -> float -> unit
  = "c_ba_ida_calc_ic_ya_ydp"
let calc_ic_ya_yd' session ?y ?y' id tout1 =
  if Id.length id <> neqs session then
    invalid_arg "calc_ic_ya_yd': variable type array (of type VarTypes.t) has incorrect length";
  let len = function
    | None -> neqs session
    | Some x -> Carray.length x
  in
  if len y <> neqs session then
    invalid_arg "calc_ic_ya_yd': ~y has incorrect length";
  if len y' <> neqs session then
    invalid_arg "calc_ic_ya_yd': ~y' has incorrect length";
  c_calc_ic_ya_yd' session y y' id tout1
