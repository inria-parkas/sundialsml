(***********************************************************************)
(*                                                                     *)
(*              Ocaml interface to Sundials CVODE solver               *)
(*                                                                     *)
(*           Timothy Bourke (INRIA) and Marc Pouzet (LIENS)            *)
(*                                                                     *)
(*  Copyright 2013 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under the terms of the GNU Library General Public License, with    *)
(*  the special exception on linking described in file LICENSE.        *)
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

type ida_mem
type ida_file
type session = {
        ida        : ida_mem;
        neqs       : int;
        nroots     : int;
        err_file   : ida_file;

        (* Temporary storage for exceptions raised within callbacks.  Not
           to be touched from the OCaml side.  *)
        mutable exn_temp   : exn option;

        mutable resfn      : float -> val_array -> der_array -> val_array
                             -> unit;
        mutable rootsfn    : float -> val_array -> der_array -> root_val_array
                             -> unit;
        mutable errh       : error_details -> unit;
        mutable errw       : val_array -> nvec -> unit;
        mutable jacfn      : triple_tmp jacobian_arg -> Densematrix.t -> unit;
        mutable bandjacfn  : triple_tmp jacobian_arg -> int -> int
                               -> Bandmatrix.t -> unit;
        mutable presetupfn : triple_tmp jacobian_arg -> unit;
        mutable presolvefn : single_tmp jacobian_arg -> val_array -> val_array
                               -> float -> unit;
        mutable jactimesfn : double_tmp jacobian_arg -> val_array -> val_array
                               -> unit;
      }

external session_finalize : session -> unit
  = "c_ida_session_finalize"

(* FIXME: isn't it better to separate out IDARootInit(), since it's
   optional in the C interface?  *)
external external_init
  : linear_solver -> nvec -> nvec -> int -> int -> float -> (ida_mem * ida_file)
  = "c_ba_ida_init_bytecode" "c_ba_ida_init"

let init' linsolv resfn (nroots, roots) y y' t0 =
  let neqs = Sundials.Carray.length y in
  (* IDA doesn't check if y and y' have the same length, and corrupt memory if
   * they aren't.  *)
  if neqs <> Sundials.Carray.length y' then
    raise (Invalid_argument "y and y' have inconsistent sizes");
  let (ida_mem, err_file) =
    external_init linsolv y y' neqs nroots t0 in
  let session = { ida        = ida_mem;
                  neqs       = neqs;
                  nroots     = nroots;
                  err_file   = err_file;
                  exn_temp   = None;
                  resfn      = resfn;
                  rootsfn    = roots;
                  errh       = (fun _ -> ());
                  errw       = (fun _ _ -> ());
                  jacfn      = (fun _ _ -> ());
                  bandjacfn  = (fun _ _ _ _ -> ());
                  presetupfn = (fun _ -> ());
                  presolvefn = (fun _ _ _ _ -> ());
                  jactimesfn = (fun _ _ _ -> ());
                }
  in
  Gc.finalise session_finalize session;
  session

let init linsolv resfn roots y yp = init' linsolv resfn roots y yp 0.
(*
  module SessionTable : sig
  val proto_init :
  lmm
  -> iter
  -> (float -> val_array -> der_array -> unit)
  -> (int * (float -> val_array -> root_val_array -> unit))
  -> nvec
  -> float
  -> session
  end = 
  struct

  let session_table = ref (Weak.create 10 : session Weak.t)

  let add_session ida neqs nroots err_file =
  let length = Weak.length !session_table in
  let rec find_next i =
  if i < length
  then (if Weak.check !session_table i then find_next (i + 1) else i)
  else
  let session_table' = Weak.create (2 * length) in
  Weak.blit !session_table 0 session_table' 0 length;
  session_table := session_table';
  i in
  let idx = find_next 0 in
  let session = {
  ida      = ida;
  neqs       = neqs;
  nroots     = nroots;
  err_file   = err_file;

  resfn      = (fun _ _ _ -> ());
  rootsfn    = (fun _ _ _ -> ());
  errh       = (fun _ -> ());
  errw       = (fun _ _ -> ());
  jacfn      = (fun _ _ -> ());
  bandjacfn  = (fun _ _ _ _ -> ());
  presetupfn = (fun _ _ _ -> false);
  presolvefn = (fun _ _ _ -> ());
  jactimesfn = (fun _ _ _ -> ());
  } in
  Weak.set !session_table idx (Some session);
  session

  let get_session idx =
  match Weak.get !session_table idx with
  None -> raise Not_found
  | Some s -> s

  let session_resfn idx      = (get_session idx).resfn
  let session_errh idx       = (get_session idx).errh
  let session_errw idx       = (get_session idx).errw
  let session_jacfn idx      = (get_session idx).jacfn
  let session_bandjacfn idx  = (get_session idx).bandjacfn
  let session_presetupfn idx = (get_session idx).presetupfn
  let session_presolvefn idx = (get_session idx).presolvefn
  let session_jactimesfn idx = (get_session idx).jactimesfn

  let _ = Callback.register "c_ba_ida_ml_session"    get_session;
  Callback.register "c_ba_ida_ml_resfn"      session_resfn;
  Callback.register "c_ba_ida_ml_errh"       session_errh;
  Callback.register "c_ba_ida_ml_errw"       session_errw;
  Callback.register "c_ba_ida_ml_jacfn"      session_jacfn;
  Callback.register "c_ba_ida_ml_bandjacfn"  session_bandjacfn;
  Callback.register "c_ba_ida_ml_presetupfn" session_presetupfn;
  Callback.register "c_ba_ida_ml_presolvefn" session_presolvefn;
  Callback.register "c_ba_ida_ml_jactimesfn" session_jactimesfn

  let proto_init lmm iter f (num_roots, roots) y0 t0 =
  let num_eqs = Sundials.Carray.length y0 in
  let ida, errfile = external_init lmm iter y0 num_eqs num_roots t0 in
  let s = add_session ida num_eqs num_roots errfile in
  Gc.finalise session_finalize s;
  s.resfn <- f;
  s.rootsfn <- roots;
  set_user_data s;
  s

  end

  let init' = SessionTable.proto_init
  let init lmm iter f roots n_y0 =
  SessionTable.proto_init lmm iter f roots n_y0 0.0
 *)

let nroots { nroots } = nroots
let neqs { neqs } = neqs

external reinit
    : session -> float -> val_array -> der_array -> unit
    = "c_ba_ida_reinit"

external sv_tolerances  : session -> float -> nvec -> unit
    = "c_ba_ida_sv_tolerances"
external ss_tolerances  : session -> float -> float -> unit
    = "c_ida_ss_tolerances"
external wf_tolerances  : session -> unit
    = "c_ba_ida_wf_tolerances"

let wf_tolerances s ferrw =
  s.errw <- ferrw;
  wf_tolerances s

external get_root_info  : session -> root_array -> unit
    = "c_ida_get_root_info"

external normal
    : session -> float -> val_array -> der_array -> float * solver_result
    = "c_ba_ida_normal"

external one_step
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
  set_root_direction' s (RootDirs.create' (nroots s) rda)

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
      s.jacfn <- (fun _ _ -> ());
      clear_dense_jac_fn s

    external set_band_jac_fn   : session -> unit
        = "c_ba_ida_dls_set_band_jac_fn"

    let set_band_jac_fn s fbandjacfn =
      s.bandjacfn <- fbandjacfn;
      set_band_jac_fn s

    external clear_band_jac_fn : session -> unit
        = "c_ba_ida_dls_clear_band_jac_fn"

    let clear_band_jac_fn s =
      s.bandjacfn <- (fun _ _ _ _ -> ());
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
    type gramschmidt_type =
      | ModifiedGS
      | ClassicalGS

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

    external set_gs_type : session -> gramschmidt_type -> unit
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
    let length = Carray.length

    let get a i = constraint_type_of_float a.{i}
    let set a i x = a.{i} <- float_of_constraint_type x
    let fill a t =
      let x = float_of_constraint_type t in
      Carray.fill a x

    let blit a b = Carray.blit a b
  end

external set_constraints : session -> Constraints.t -> unit
  = "c_ba_ida_set_constraints"

module Id =
  struct
    type t = nvec
    type component_type = Algebraic | Differential

    let component_type_of_float = function
      | 0.0 -> Algebraic
      | 1.0 -> Differential
      | f -> raise (Invalid_argument
                      ("invalid component type: " ^ string_of_float f))
    let float_of_component_type = function
      | Algebraic -> 0.0
      | Differential -> 1.0

    let create = Carray.create
    let length = Carray.length

    let get a i = component_type_of_float a.{i}
    let set a i x = a.{i} <- float_of_component_type x
    let set_algebraic a i = set a i Algebraic
    and set_differential a i = set a i Differential
    let fill a t =
      let x = float_of_component_type t in
      Carray.fill a x

    let blit a b = Carray.blit a b
  end

external calc_ic_y : session -> float -> unit
  = "c_ba_ida_calc_ic_y"
external c_calc_ic_ya_yd' : session -> Id.t -> float -> unit
  = "c_ba_ida_calc_ic_ya_ydp"
let calc_ic_ya_yd' session id tout1 =
  if Id.length id <> neqs session then
    raise (Invalid_argument ("length of component type array does not match number of equations"));
  c_calc_ic_ya_yd' session id tout1
