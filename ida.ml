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
include Ida_impl

(*
 * NB: The order of variant constructors and record fields is important!
 *     If these types are changed or augmented, the corresponding declarations
 *     in ida_ml.h (and code in ida_ml.c) must also be updated.
 *)

(* Solver exceptions *)
exception IllInput
exception TooMuchWork
exception TooMuchAccuracy
exception ErrFailure
exception ConvergenceFailure
exception LinearInitFailure
exception LinearSetupFailure
exception LinearSolveFailure
exception ResFuncFailure
exception FirstResFuncFailure
exception RepeatedResFuncFailure
exception RootFuncFailure
exception ConstraintFailure

(* Initial condition calculator exceptions *)
exception LinesearchFailure
exception NoRecovery
exception BadEwt

(* get_dky exceptions *)
exception BadK
exception BadT

(* SetId exceptions *)
exception IdNotSet

let no_roots = (0, dummy_rootsfn)

type integrator_stats = {
    num_steps : int;
    num_res_evals : int;
    num_lin_solv_setups : int;
    num_err_test_fails : int;
    last_order : int;
    current_order : int;
    actual_init_step : float;
    last_step : float;
    current_step : float;
    current_time : float
  }

module VarId =
  struct
    let algebraic = 0.0
    let differential = 1.0

    type t = Algebraic | Differential

    let of_float x =
      if x = algebraic then Algebraic
      else if x = differential then Differential
      else invalid_arg ("invalid component type: " ^ string_of_float x)

    let to_float = function
      | Algebraic -> algebraic
      | Differential -> differential
  end

external c_root_init : ('a, 'k) session -> int -> unit
    = "c_ida_root_init"

let root_init session (nroots, rootsfn) =
  c_root_init session nroots;
  session.rootsfn <- rootsfn

module Dls =
  struct
    include DlsTypes

    external c_dls_dense : serial_session -> int -> bool -> unit
      = "c_ida_dls_dense"

    external c_dls_lapack_dense : serial_session -> int -> bool -> unit
      = "c_ida_dls_lapack_dense"

    external c_dls_band : serial_session -> int -> int -> int -> bool -> unit
      = "c_ida_dls_band"

    external c_dls_lapack_band : serial_session -> int -> int -> int -> bool
                               -> unit
      = "c_ida_dls_lapack_band"

    external set_dense_jac_fn : serial_session -> unit
        = "c_ida_dls_set_dense_jac_fn"

    let dense ?jac () session nv nv' =
      let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
      (session.ls_callbacks <-
        match jac with
        | None -> NoCallbacks
        | Some f -> DenseCallback { jacfn = f; dmat = None });
      c_dls_dense session neqs (jac <> None)

    let lapack_dense ?jac () session nv nv' =
      let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
      (session.ls_callbacks <-
        match jac with
        | None -> NoCallbacks
        | Some f -> DenseCallback { jacfn = f; dmat = None });
      c_dls_lapack_dense session neqs (jac <> None)

    let band ?jac p session nv nv' =
      let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
      (session.ls_callbacks <-
        match jac with
        | None -> NoCallbacks
        | Some f -> BandCallback { bjacfn = f; bmat = None });
      c_dls_band session neqs p.mupper p.mlower (jac <> None)

    let lapack_band ?jac p session nv nv' =
      let neqs = Sundials.RealArray.length (Nvector.unwrap nv) in
      (session.ls_callbacks <-
        match jac with
        | None -> NoCallbacks
        | Some f -> BandCallback { bjacfn = f; bmat = None });
      c_dls_lapack_band session neqs p.mupper p.mlower (jac <> None)


    let invalidate_callback (type d) (type k) (session : (d, k) session) =
      match session.ls_callbacks with
      | DenseCallback ({ dmat = Some d } as cb) ->
          Dls.DenseMatrix.invalidate d;
          cb.dmat <- None
      | BandCallback  ({ bmat = Some d } as cb) ->
          Dls.BandMatrix.invalidate d;
          cb.bmat <- None
      | _ -> ()

    let set_dense_jac_fn s fjacfn =
      invalidate_callback s;
      s.ls_callbacks <- DenseCallback { jacfn = fjacfn; dmat = None };
      set_dense_jac_fn s

    external clear_dense_jac_fn : serial_session -> unit
        = "c_ida_dls_clear_dense_jac_fn"

    let clear_dense_jac_fn s =
      match s.ls_callbacks with
      | DenseCallback _ -> (invalidate_callback s;
                            s.ls_callbacks <- NoCallbacks;
                            clear_dense_jac_fn s)
      | _ -> failwith "dense linear solver not in use"

    external set_band_jac_fn : serial_session -> unit
        = "c_ida_dls_set_band_jac_fn"

    let set_band_jac_fn s fbandjacfn =
      invalidate_callback s;
      s.ls_callbacks <- BandCallback { bjacfn = fbandjacfn; bmat = None };
      set_band_jac_fn s

    external clear_band_jac_fn : serial_session -> unit
        = "c_ida_dls_clear_band_jac_fn"

    let clear_band_jac_fn s =
      match s.ls_callbacks with
      | BandCallback _ -> (invalidate_callback s;
                           s.ls_callbacks <- NoCallbacks;
                           clear_band_jac_fn s)
      | _ -> failwith "banded linear solver not in use"

    external get_work_space : serial_session -> int * int
        = "c_ida_dls_get_work_space"

    external get_num_jac_evals : serial_session -> int
        = "c_ida_dls_get_num_jac_evals"

    external get_num_res_evals : serial_session -> int
        = "c_ida_dls_get_num_res_evals"
  end

module Spils =
  struct
    include SpilsTypes

    external c_spgmr
      : ('a, 'k) session -> int -> unit
      = "c_ida_spils_spgmr"

    external c_spbcg
      : ('a, 'k) session -> int -> unit
      = "c_ida_spils_spbcg"

    external c_sptfqmr
      : ('a, 'k) session -> int -> unit
      = "c_ida_spils_sptfqmr"

    external c_set_preconditioner
      : ('a, 'k) session -> bool -> unit
      = "c_ida_spils_set_preconditioner"

    external c_set_jac_times_vec_fn : ('a, 'k) session -> bool -> unit
      = "c_ida_spils_set_jac_times_vec_fn"

    external c_set_max_restarts : ('a, 'k) session -> int -> unit
      = "c_ida_spils_set_max_restarts"

    let init_preconditioner solve setup jac_times session nv nv' =
      c_set_preconditioner session (setup <> None);
      c_set_jac_times_vec_fn session (jac_times <> None);
      session.ls_callbacks <-
        SpilsCallback { prec_solve_fn = solve;
                        prec_setup_fn = setup;
                        jac_times_vec_fn = jac_times }

    let prec_none = InternalPrecNone
    let prec_left ?setup ?jac_times_vec solve =
      InternalPrecLeft (init_preconditioner solve setup jac_times_vec)

    let init_spils init maxl prec session nv nv' =
      match prec with
      | InternalPrecNone -> init session maxl
      | InternalPrecLeft set_prec ->
        init session maxl;
        set_prec session nv nv'

    let spgmr ?(maxl=0) ?max_restarts prec session nv nv' =
      init_spils c_spgmr maxl prec session nv nv';
      (match max_restarts with
       | Some m -> c_set_max_restarts session m
       | None -> ())

    let spbcg ?(maxl=0) prec session nv nv' =
      init_spils c_spbcg maxl prec session nv nv'

    let sptfqmr ?(maxl=0) prec session nv nv' =
      init_spils c_sptfqmr maxl prec session nv nv'

    let set_preconditioner s ?setup solve =
      match s.ls_callbacks with
      | SpilsCallback cbs ->
        c_set_preconditioner s (setup <> None);
        s.ls_callbacks <- SpilsCallback { cbs with
                                          prec_setup_fn = setup;
                                          prec_solve_fn = solve }
      | _ -> failwith "spils solver not in use"

    let set_jac_times_vec_fn s f =
      match s.ls_callbacks with
      | SpilsCallback cbs ->
        c_set_jac_times_vec_fn s true;
        s.ls_callbacks <- SpilsCallback { cbs with jac_times_vec_fn = Some f }
      | _ -> failwith "spils solver not in use"

    let clear_jac_times_vec_fn s =
      match s.ls_callbacks with
      | SpilsCallback cbs ->
        c_set_jac_times_vec_fn s false;
        s.ls_callbacks <- SpilsCallback { cbs with jac_times_vec_fn = None }
      | _ -> failwith "spils solver not in use"

    external set_gs_type : ('a, 'k) session -> Spils.gramschmidt_type -> unit
        = "c_ida_spils_set_gs_type"

    external set_eps_lin            : ('a, 'k) session -> float -> unit
        = "c_ida_spils_set_eps_lin"

    external set_maxl               : ('a, 'k) session -> int -> unit
        = "c_ida_spils_set_maxl"

    external get_num_lin_iters      : ('a, 'k) session -> int
        = "c_ida_spils_get_num_lin_iters"

    external get_num_conv_fails     : ('a, 'k) session -> int
        = "c_ida_spils_get_num_conv_fails"

    external get_work_space         : ('a, 'k) session -> int * int
        = "c_ida_spils_get_work_space"

    external get_num_prec_evals     : ('a, 'k) session -> int
        = "c_ida_spils_get_num_prec_evals"

    external get_num_prec_solves    : ('a, 'k) session -> int
        = "c_ida_spils_get_num_prec_solves"

    external get_num_jtimes_evals   : ('a, 'k) session -> int
        = "c_ida_spils_get_num_jtimes_evals"

    external get_num_res_evals      : ('a, 'k) session -> int
        = "c_ida_spils_get_num_res_evals"

  end

module Alternate =
  struct
    include AlternateTypes

    external c_set_alternate
      : ('data, 'kind) session -> bool -> bool -> unit
      = "c_ida_set_alternate"

    let make f s nv nv' =
      let { linit; lsetup; lsolve } as cb = f s nv nv' in
      c_set_alternate s (linit <> None) (lsetup <> None);
      s.ls_callbacks <- AlternateCallback cb

    external get_cj : ('data, 'kind) session -> float = "c_ida_get_cj"
    external get_cjratio : ('data, 'kind) session -> float = "c_ida_get_cjratio"
  end


let set_linear_solver session solver nv nv' =
  session.ls_callbacks <- NoCallbacks;
  solver session nv nv'

external sv_tolerances
    : ('a, 'k) session -> float -> ('a, 'k) Nvector.t -> unit
    = "c_ida_sv_tolerances"
external ss_tolerances  : ('a, 'k) session -> float -> float -> unit
    = "c_ida_ss_tolerances"
external wf_tolerances  : ('a, 'k) session -> unit
    = "c_ida_wf_tolerances"

type 'a error_fun = 'a -> 'a -> unit

type ('a, 'k) tolerance =
  | SStolerances of float * float
  | SVtolerances of float * ('a, 'k) Nvector.t
  | WFtolerances of 'a error_fun

let default_tolerances = SStolerances (1.0e-4, 1.0e-8)

let set_tolerances s tol =
  match tol with
  | SStolerances (rel, abs) -> ss_tolerances s rel abs
  | SVtolerances (rel, abs) -> (s.checkvec abs; sv_tolerances s rel abs)
  | WFtolerances ferrw -> (s.errw <- ferrw; wf_tolerances s)

external c_set_id : ('a,'k) session -> ('a,'k) Nvector.t -> unit
  = "c_ida_set_id"

let set_id s id =
  s.checkvec id;
  c_set_id s id;
  s.id_set <- true

external c_session_finalize : ('a, 'kind) session -> unit
    = "c_ida_session_finalize"

let session_finalize s =
  Dls.invalidate_callback s;
  c_session_finalize s

external c_init : ('a, 'k) session Weak.t -> float
                  -> ('a, 'k) Nvector.t -> ('a, 'k) Nvector.t
                  -> (ida_mem * c_weak_ref * ida_file)
    = "c_ida_init"

let init linsolv tol resfn ?varid ?(roots=no_roots) t0 y y' =
  let (nroots, rootsfn) = roots in
  let checkvec = Nvector.check y in
  checkvec y';
  if nroots < 0 then
    raise (Invalid_argument "number of root functions is negative");
  (* FIXME: can we check y and y' have the same length, at least for
     some nvector types?  *)
  let weakref = Weak.create 1 in
  let ida_mem, backref, err_file = c_init weakref t0 y y' in
  (* ida_mem and backref have to be immediately captured in a session and
     associated with the finalizer before we do anything else.  *)
  let session = { ida        = ida_mem;
                  backref    = backref;
                  nroots     = nroots;
                  err_file   = err_file;
                  checkvec   = checkvec;
                  exn_temp   = None;
                  id_set     = false;
                  resfn      = resfn;
                  rootsfn    = rootsfn;
                  errh       = dummy_errh;
                  errw       = dummy_errw;
                  ls_callbacks = NoCallbacks;
                  sensext    = NoSensExt;
                }
  in
  Gc.finalise session_finalize session;
  Weak.set weakref 0 (Some session);
  (* Now the session is safe to use.  If any of the following fails and raises
     an exception, the GC will take care of freeing ida_mem and backref.  *)
  (match varid with
     None -> ()
   | Some x -> set_id session x);
  if nroots > 0 then
    c_root_init session nroots;
  set_linear_solver session linsolv y y';
  set_tolerances session tol;
  session

let get_num_roots { nroots } = nroots

external c_reinit
    : ('a, 'k) session -> float -> ('a, 'k) Nvector.t
      -> ('a, 'k) Nvector.t -> unit
    = "c_ida_reinit"

let reinit session ?linsolv ?roots t0 y0 y'0 =
  session.checkvec y0;
  session.checkvec y'0;
  Dls.invalidate_callback session;
  c_reinit session t0 y0 y'0;
  (match linsolv with
   | None -> ()
   | Some linsolv -> set_linear_solver session linsolv y0 y'0);
  (match roots with
   | None -> ()
   | Some roots -> root_init session roots)

external get_root_info  : ('a, 'k) session -> Sundials.Roots.t -> unit
    = "c_ida_get_root_info"

type solver_result =
  | Success             (** IDA_SUCCESS *)
  | RootsFound          (** IDA_ROOT_RETURN *)
  | StopTimeReached     (** IDA_TSTOP_RETURN *)

external c_solve_normal : ('a, 'k) session -> float
                          -> ('a, 'k) Nvector.t -> ('a,'k) Nvector.t
                          -> float * solver_result
    = "c_ida_solve_normal"

let solve_normal s t y yp =
  s.checkvec y;
  s.checkvec yp;
  c_solve_normal s t y yp

external c_solve_one_step : ('a, 'k) session -> float
                            -> ('a, 'k) Nvector.t-> ('a, 'k) Nvector.t
                            -> float * solver_result
    = "c_ida_solve_one_step"

let solve_one_step s t y yp =
  s.checkvec y;
  s.checkvec yp;
  c_solve_one_step s t y yp

external c_get_dky
    : ('a, 'k) session -> float -> int -> ('a, 'k) Nvector.t -> unit
    = "c_ida_get_dky"

let get_dky s t k y =
  s.checkvec y;
  c_get_dky s t k y

external get_integrator_stats : ('a, 'k) session -> integrator_stats
    = "c_ida_get_integrator_stats"

external get_work_space         : ('a, 'k) session -> int * int
    = "c_ida_get_work_space"

external get_num_steps          : ('a, 'k) session -> int
    = "c_ida_get_num_steps"

external get_num_res_evals      : ('a, 'k) session -> int
    = "c_ida_get_num_res_evals"

external get_num_lin_solv_setups : ('a, 'k) session -> int
    = "c_ida_get_num_lin_solv_setups"

external get_num_err_test_fails : ('a, 'k) session -> int
    = "c_ida_get_num_err_test_fails"

external get_last_order         : ('a, 'k) session -> int
    = "c_ida_get_last_order"

external get_current_order      : ('a, 'k) session -> int
    = "c_ida_get_current_order"

external get_actual_init_step   : ('a, 'k) session -> float
    = "c_ida_get_actual_init_step"

external get_last_step          : ('a, 'k) session -> float
    = "c_ida_get_last_step"

external get_current_step       : ('a, 'k) session -> float
    = "c_ida_get_current_step"

external get_current_time       : ('a, 'k) session -> float
    = "c_ida_get_current_time"

let print_integrator_stats s oc =
  let stats = get_integrator_stats s
  in
    Printf.fprintf oc "num_steps = %d\n"           stats.num_steps;
    Printf.fprintf oc "num_res_evals = %d\n"       stats.num_res_evals;
    Printf.fprintf oc "num_lin_solv_setups = %d\n" stats.num_lin_solv_setups;
    Printf.fprintf oc "num_err_test_fails = %d\n"  stats.num_err_test_fails;
    Printf.fprintf oc "last_order = %d\n"          stats.last_order;
    Printf.fprintf oc "current_order = %d\n"       stats.current_order;
    Printf.fprintf oc "actual_init_step = %e\n"    stats.actual_init_step;
    Printf.fprintf oc "last_step = %e\n"           stats.last_step;
    Printf.fprintf oc "current_step = %e\n"        stats.current_step;
    Printf.fprintf oc "current_time = %e\n"        stats.current_time;

external set_error_file : ('a, 'k) session -> string -> bool -> unit
    = "c_ida_set_error_file"

external set_err_handler_fn  : ('a, 'k) session -> unit
    = "c_ida_set_err_handler_fn"

let set_err_handler_fn s ferrh =
  s.errh <- ferrh;
  set_err_handler_fn s

external clear_err_handler_fn  : ('a, 'k) session -> unit
    = "c_ida_clear_err_handler_fn"

let clear_err_handler_fn s =
  s.errh <- dummy_errh;
  clear_err_handler_fn s

external set_max_ord            : ('a, 'k) session -> int -> unit
    = "c_ida_set_max_ord"
external set_max_num_steps      : ('a, 'k) session -> int -> unit
    = "c_ida_set_max_num_steps"
external set_init_step          : ('a, 'k) session -> float -> unit
    = "c_ida_set_init_step"
external set_max_step           : ('a, 'k) session -> float -> unit
    = "c_ida_set_max_step"
external set_stop_time          : ('a, 'k) session -> float -> unit
    = "c_ida_set_stop_time"
external set_max_err_test_fails : ('a, 'k) session -> int -> unit
    = "c_ida_set_max_err_test_fails"
external set_max_nonlin_iters   : ('a, 'k) session -> int -> unit
    = "c_ida_set_max_nonlin_iters"
external set_max_conv_fails     : ('a, 'k) session -> int -> unit
    = "c_ida_set_max_conv_fails"
external set_nonlin_conv_coef   : ('a, 'k) session -> float -> unit
    = "c_ida_set_nonlin_conv_coef"

external set_root_direction'   : ('a, 'k) session -> Sundials.RootDirs.t -> unit
    = "c_ida_set_root_direction"

let set_root_direction s rda =
  set_root_direction' s (Sundials.RootDirs.copy (get_num_roots s) rda)

let set_all_root_directions s rd =
  set_root_direction' s (Sundials.RootDirs.make (get_num_roots s) rd)

external set_no_inactive_root_warn      : ('a, 'k) session -> unit
    = "c_ida_set_no_inactive_root_warn"
(*
   IDAGetNumStabLimOrderReds appears in the sundials 2.5.0 manual on
   p.52 but there's no such function in the implementation.  It's
   probably a leftover from earlier versions or something.

external get_num_stab_lim_order_reds    : ('a, 'k) session -> int
    = "c_ida_get_num_stab_lim_order_reds"
*)
external get_tol_scale_factor           : ('a, 'k) session -> float
    = "c_ida_get_tol_scale_factor"

external c_get_err_weights : ('a, 'k) session -> ('a, 'k) Nvector.t -> unit
    = "c_ida_get_err_weights"

let get_err_weights s ew =
  s.checkvec ew;
  c_get_err_weights s ew

external c_get_est_local_errors : ('a, 'k) session -> ('a, 'k) Nvector.t -> unit
    = "c_ida_get_est_local_errors"

let get_est_local_errors s ew =
  s.checkvec ew;
  c_get_est_local_errors s ew

external get_num_nonlin_solv_iters      : ('a, 'k) session -> int
    = "c_ida_get_num_nonlin_solv_iters"

external get_num_nonlin_solv_conv_fails : ('a, 'k) session -> int
    = "c_ida_get_num_nonlin_solv_conv_fails"

external get_nonlin_solv_stats          : ('a, 'k) session -> int * int
    = "c_ida_get_nonlin_solv_stats"

external get_num_g_evals                : ('a, 'k) session -> int
    = "c_ida_get_num_g_evals"

external c_set_constraints : ('a,'k) session -> ('a,'k) Nvector.t -> unit
  = "c_ida_set_constraints"

let set_constraints s nv =
  s.checkvec nv;
  c_set_constraints s nv

external c_set_suppress_alg : ('a,'k) session -> bool -> unit
  = "c_ida_set_suppress_alg"

let set_suppress_alg s ?varid v =
  (match varid with
   | None -> if v && not s.id_set then raise IdNotSet
   | Some x -> set_id s x);
  c_set_suppress_alg s v

external get_num_backtrack_ops : ('a,'k) session -> int
  = "c_ida_get_num_backtrack_ops"

external c_calc_ic_y : ('a,'k) session -> ('a,'k) Nvector.t option
                       -> float -> unit
  = "c_ida_calc_ic_y"

let calc_ic_y session ?y tout1 =
  (match y with None -> () | Some x -> session.checkvec x);
  c_calc_ic_y session y tout1

external c_calc_ic_ya_yd' :
  ('a,'k) session -> ('a,'k) Nvector.t option -> ('a,'k) Nvector.t option
  -> float -> unit
  = "c_ida_calc_ic_ya_ydp"

let calc_ic_ya_yd' session ?y ?y' ?varid tout1 =
  (match y with None -> () | Some x -> session.checkvec x);
  (match y' with None -> () | Some x -> session.checkvec x);
  (match varid with
   | None -> if not session.id_set then raise IdNotSet
   | Some x -> set_id session x);
  c_calc_ic_ya_yd' session y y' tout1

(* Callbacks *)

let call_errw session y ewt =
  let session = read_weak_ref session in
  try session.errw y ewt; 0
  with
  | Sundials.NonPositiveEwt -> -1
  | e -> (session.exn_temp <- Some e; -1)

let call_errh session details =
  let session = read_weak_ref session in
  try session.errh details
  with e ->
    prerr_endline ("Warning: error handler function raised an exception.  " ^
                   "This exception will not be propagated: " ^
                   Printexc.to_string e)

(* Let C code know about some of the values in this module.  *)
external c_init_module : 'fcns -> exn array -> unit =
  "c_ida_init_module"

let _ =
  c_init_module
    (* Functions must be listed in the same order as
       callback_index in ida_ml.c.  *)
    (call_errh, call_errw)

    (* Exceptions must be listed in the same order as
       ida_exn_index.  *)
    [|IllInput;
      TooMuchWork;
      TooMuchAccuracy;
      ErrFailure;
      ConvergenceFailure;
      LinearInitFailure;
      LinearSetupFailure;
      LinearSolveFailure;
      ResFuncFailure;
      FirstResFuncFailure;
      RepeatedResFuncFailure;
      RootFuncFailure;
      ConstraintFailure;

      LinesearchFailure;
      NoRecovery;
      BadEwt;

      BadK;
      BadT;
    |]
