(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2020 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

open Sundials

include Sundials_NonlinearSolver_impl

type ('data, 'kind, 's, 'v) t
    = ('data, 'kind, 's, 'v) Sundials_NonlinearSolver_impl.nonlinear_solver

let check_compat () =
  if Sundials_impl.Version.in_compat_mode2_3
  then raise Config.NotImplementedBySundialsVersion

let weak_wrap x =
  let wx = Weak.create 1 in
  Weak.set wx 0 (Some x);
  wx

external c_init : ('d, 'k, 's, 'v) cptr -> unit
  = "sunml_nlsolver_init"

external c_setup :
  ('d, 'k, 's, [`Nvec]) cptr
  -> ('d, 'k) Nvector.t
  -> 's
  -> unit
  = "sunml_nlsolver_setup"

external c_setup_sens :
  ('d, 'k, 's, [`Sens]) cptr
  -> ('d, 'k) Senswrapper.t
  -> 's
  -> unit
  = "sunml_nlsolver_setup_sens"

external c_solve :
  ('d, 'k, 's, [`Nvec]) cptr
  -> (('d, 'k) Nvector.t * ('d, 'k) Nvector.t * ('d, 'k) Nvector.t * float * bool)
  -> 's
  -> unit
  = "sunml_nlsolver_solve"

external c_solve_sens :
  ('d, 'k, 's, [`Sens]) cptr
  -> (('d, 'k) Senswrapper.t * ('d, 'k) Senswrapper.t * ('d, 'k) Senswrapper.t * float * bool)
  -> 's
  -> unit
  = "sunml_nlsolver_solve_sens"

external c_set_sys_fn      : ('d, 'k, 's, [`Nvec]) cptr -> unit
  = "sunml_nlsolver_set_sys_fn"

external c_set_sys_fn_sens : ('d, 'k, 's, [`Sens]) cptr -> unit
  = "sunml_nlsolver_set_sys_fn_sens"

external c_set_lsetup_fn   : ('d, 'k, 's, 'v) cptr -> unit
  = "sunml_nlsolver_set_lsetup_fn"

external c_set_lsolve_fn   : ('d, 'k, 's, [`Nvec]) cptr -> unit
  = "sunml_nlsolver_set_lsolve_fn"

external c_set_lsolve_fn_sens : ('d, 'k, 's, [`Sens]) cptr -> unit
  = "sunml_nlsolver_set_lsolve_fn_sens"

external c_set_convtest_fn
  :    ('d, 'k, 's, [`Nvec]) cptr
    -> 's convtest_callback Sundials.cfun
    -> unit
  = "sunml_nlsolver_set_convtest_fn"

external c_set_convtest_fn_callback
  : ('d, 'k, 's, [`Nvec]) cptr -> unit
  = "sunml_nlsolver_set_convtest_fn_callback"

external c_set_convtest_fn_sens
  :    ('d, 'k, 's, [`Sens]) cptr
    -> 's convtest_callback_sens Sundials.cfun
    -> unit
  = "sunml_nlsolver_set_convtest_fn" (* Same as c_set_convtest_fn *)

external c_set_convtest_fn_callback_sens
  : ('d, 'k, 's, [`Sens]) cptr -> unit
  = "sunml_nlsolver_set_convtest_fn_callback_sens"

external c_set_max_iters : ('d, 'k, 's, 'v) cptr -> int -> unit
  = "sunml_nlsolver_set_max_iters"

external c_get_num_iters : ('d, 'k, 's, 'v) cptr -> int
  = "sunml_nlsolver_get_num_iters"

external c_get_cur_iter  : ('d, 'k, 's, 'v) cptr -> int
  = "sunml_nlsolver_get_cur_iter"

external c_get_num_conv_fails : ('d, 'k, 's, 'v) cptr -> int
  = "sunml_nlsolver_get_num_conv_fails"

external c_set_info_file_fixedpoint : ('d, 'k, 's, 'v) cptr -> Logfile.t -> unit
  = "sunml_nlsolver_set_info_file_fixedpoint"

external c_set_info_file_newton : ('d, 'k, 's, 'v) cptr -> Logfile.t -> unit
  = "sunml_nlsolver_set_info_file_newton"

external c_set_print_level_fixedpoint : ('d, 'k, 's, 'v) cptr -> int -> unit
  = "sunml_nlsolver_set_print_level_fixedpoint"

external c_set_print_level_newton : ('d, 'k, 's, 'v) cptr -> int -> unit
  = "sunml_nlsolver_set_print_level_newton"

(* - - - OCaml invoking init/setup/solve - - - *)

let uw = Nvector.unwrap

let init (type v) ({ rawptr; solver } : ('d, 'k, 's, v) t) =
  check_compat ();
  match solver with
  | CustomSolver     (_, { init = Some f }) -> f ()          (* O/Onls *)
  | CustomSolverSens (_, { init = Some f }) -> f ()          (* O/Onls *)
  | CustomSolver _ -> ()
  | CustomSolverSens _   -> ()
  | FixedPointSolver _ -> c_init rawptr                      (* O/Cnls *)
  | NewtonSolver _ -> c_init rawptr
  | FixedPointSolverSens _ -> c_init rawptr
  | NewtonSolverSens _ -> c_init rawptr

let setup { rawptr; solver } ~y s =
  check_compat ();
  match solver with
  | CustomSolver     (_, { setup = Some f }) -> f (uw y) s   (* O/Onls *)
  | CustomSolver     (_, { setup = None   }) -> ()
  | FixedPointSolver _ | NewtonSolver _
      -> c_setup rawptr y s                                  (* O/Cnls *)

let solve { rawptr; solver } ~y0 ~ycor ~w tol callLSetup s =
  check_compat ();
  match solver with
  | CustomSolver (_, { solve = f })                          (* O/Onls *)
      -> f (uw y0) (uw ycor) (uw w) tol callLSetup s
  | FixedPointSolver _ | NewtonSolver _
      -> c_solve rawptr (y0, ycor, w, tol, callLSetup) s     (* O/Cnls *)

(* - - - OCaml callback configuration - - - *)

let set_sys_fn { rawptr; solver; _ } cbf =
  check_compat ();
  match solver with
  | CustomSolver (_, { set_sys_fn = set })                   (* O/Onls *)
      -> set (fun y fg -> cbf (uw y) (uw fg))
  | FixedPointSolver (callbacks, _) | NewtonSolver callbacks (* O/Cnls *)
      -> callbacks.sysfn <- cbf;
         c_set_sys_fn rawptr

let set_lsetup_fn (type v) ({ rawptr; solver; _ } : ('d, 'k, 's, v) t) cbf =
  check_compat ();
  match solver with
  | CustomSolver (_, { set_lsetup_fn = Some set }) -> set cbf (* O/Onls *)
  | CustomSolver (_, { set_lsetup_fn = None }) -> ()
  | CustomSolverSens (_, { set_lsetup_fn = Some set }) -> set cbf
  | CustomSolverSens (_, { set_lsetup_fn = None }) -> ()
  | FixedPointSolver (callbacks, _) ->                        (* O/Cnls *)
      callbacks.lsetupfn <- cbf;
      c_set_lsetup_fn rawptr
  | NewtonSolver callbacks ->
      callbacks.lsetupfn <- cbf;
      c_set_lsetup_fn rawptr
  | FixedPointSolverSens (callbacks, _) ->
      callbacks.lsetupfn <- cbf;
      c_set_lsetup_fn rawptr
  | NewtonSolverSens callbacks ->
      callbacks.lsetupfn <- cbf;
      c_set_lsetup_fn rawptr

let set_lsolve_fn { rawptr; solver; _ } cbf =
  check_compat ();
  match solver with
  | CustomSolver (_, { set_lsolve_fn = Some set })            (* O/Onls *)
      -> set (fun b -> cbf (uw b))
  | CustomSolver (_, { set_lsolve_fn = None }) -> ()
  | FixedPointSolver (callbacks, _) | NewtonSolver callbacks  (* O/Cnls *)
      -> callbacks.lsolvefn <- cbf;
         c_set_lsolve_fn rawptr

let assert_not_oconvtestfn (type d s) (ctfn : (d, s, [`Nvec]) convtestfn) =
  match ctfn with
  | CConvTest cfn -> CConvTest cfn
  | OConvTest _ -> invalid_arg "conftestfn is oconvtestfn"

let set_convtest_fn (type d k s)
                    ({ rawptr; solver; _ } : (d, k, s, [`Nvec]) t)
                    (ctfn : (d, s, [`Nvec]) convtestfn) =
  check_compat ();
  match solver with
  | CustomSolver (_, { set_convtest_fn = Some set }) -> set ctfn  (* O/Onls *)
  | CustomSolver _ -> ()
  | FixedPointSolver (callbacks, _) | NewtonSolver callbacks      (* O/Cnls *)
      -> (match ctfn with
        | CConvTest cfun ->
              callbacks.convtestfn <- empty_convtestfn;
              c_set_convtest_fn rawptr cfun
          | OConvTest ofn ->
              callbacks.convtestfn <- ofn;
              c_set_convtest_fn_callback rawptr)

module Sens = struct (* {{{ *)

  let setup { rawptr; solver } ~y s =
    check_compat ();
    match solver with
    | CustomSolverSens (_, { setup = Some f }) -> f y s        (* O/Onls *)
    | CustomSolverSens (_, { setup = None   }) -> ()
    | FixedPointSolverSens _ | NewtonSolverSens _
        -> c_setup_sens rawptr y s                             (* O/Cnls *)

  let solve { rawptr; solver } ~y0 ~ycor ~w tol callLSetup s =
    check_compat ();
    match solver with
    | CustomSolverSens (_, { solve = f })                       (* O/Onls *)
        -> f y0 ycor w tol callLSetup s
    | FixedPointSolverSens _ | NewtonSolverSens _
        -> c_solve_sens rawptr (y0, ycor, w, tol, callLSetup) s (* O/Cnls *)

  let set_sys_fn { rawptr; solver; _ } cbf =
    check_compat ();
    match solver with
    | CustomSolverSens (_, { set_sys_fn = set })                (* O/Onls *)
        -> set cbf
    | FixedPointSolverSens (callbacks, _) | NewtonSolverSens callbacks
        -> callbacks.sysfn <- cbf;
           c_set_sys_fn_sens rawptr

  let set_lsolve_fn { rawptr; solver; _ } cbf =
    check_compat ();
    match solver with
    | CustomSolverSens (_, { set_lsolve_fn = Some set })        (* O/Onls *)
        -> set (fun b -> cbf b)
    | CustomSolverSens (_, { set_lsolve_fn = None }) -> ()
    | FixedPointSolverSens (callbacks, _) | NewtonSolverSens callbacks
        -> callbacks.lsolvefn <- cbf;
           c_set_lsolve_fn_sens rawptr

  let assert_not_oconvtestfn (type d s) (ctfn : (d, s, [`Sens]) convtestfn) =
    match ctfn with
    | CSensConvTest cfn -> CSensConvTest cfn
    | OConvTest _ -> invalid_arg "conftestfn is oconvtestfn"

  let set_convtest_fn (type d k s)
                      ({ rawptr; solver; _ } : (d, k, s, [`Sens]) t)
                      (ctfn : ((d, k) Senswrapper.t, s, [`Sens]) convtestfn) =
    check_compat ();
    match solver with
    | CustomSolverSens (_, { set_convtest_fn = Some set }) -> set ctfn
    | CustomSolverSens _ -> ()
    | FixedPointSolverSens (callbacks, _) | NewtonSolverSens callbacks  (* O/Cnls *)
        -> (match ctfn with
            | CSensConvTest cfun ->
                callbacks.convtestfn <- empty_convtestfn;
                c_set_convtest_fn_sens rawptr cfun
            | OConvTest ofn ->
                callbacks.convtestfn <- ofn;
                c_set_convtest_fn_callback_sens rawptr)

end (* }}} *)

let set_max_iters (type d k s v) ({ rawptr; solver; _ } : (d, k, s, v) t) i =
  check_compat ();
  match solver with
  | CustomSolver     (_, { set_max_iters = Some f }) -> f i
  | CustomSolverSens (_, { set_max_iters = Some f }) -> f i
  | CustomSolver _ -> ()
  | CustomSolverSens _ -> ()
  | FixedPointSolverSens _ -> c_set_max_iters rawptr i
  | NewtonSolverSens _ -> c_set_max_iters rawptr i
  | FixedPointSolver _ -> c_set_max_iters rawptr i
  | NewtonSolver _ -> c_set_max_iters rawptr i

let set_print_level (type d k s v) ({ rawptr; solver; _ } : (d, k, s, v) t) level =
  if Sundials_impl.Version.lt530
    then raise Config.NotImplementedBySundialsVersion;
  let level = if level then 1 else 0 in
  match solver with
  | CustomSolver     (_, { set_print_level = Some f }) -> f level
  | CustomSolverSens (_, { set_print_level = Some f }) -> f level
  | CustomSolver _ -> ()
  | CustomSolverSens _ -> ()
  | FixedPointSolverSens _ -> c_set_print_level_fixedpoint rawptr level
  | NewtonSolverSens _ -> c_set_print_level_newton rawptr level
  | FixedPointSolver _ -> c_set_print_level_fixedpoint rawptr level
  | NewtonSolver _ -> c_set_print_level_newton rawptr level

let set_info_file (type d k s v)
                  ({ rawptr; solver; _ } as s : (d, k, s, v) t) ?print_level file =
  if Sundials_impl.Version.lt530
    then raise Config.NotImplementedBySundialsVersion;
  s.info_file <- Some file;
  (match solver with
   | CustomSolver     (_, { set_info_file = Some f }) -> f file
   | CustomSolverSens (_, { set_info_file = Some f }) -> f file
   | CustomSolver _ -> ()
   | CustomSolverSens _ -> ()
   | FixedPointSolverSens _ -> c_set_info_file_fixedpoint rawptr file
   | FixedPointSolver _ -> c_set_info_file_fixedpoint rawptr file
   | NewtonSolverSens _ -> c_set_info_file_newton rawptr file
   | NewtonSolver _ -> c_set_info_file_newton rawptr file);
  (match print_level with None -> () | Some level -> set_print_level s level)

let get_num_iters (type d k s v) ({ rawptr; solver; _ } : (d, k, s, v) t) =
  check_compat ();
  match solver with
  | CustomSolver     (_, { get_num_iters = Some f }) -> f ()
  | CustomSolverSens (_, { get_num_iters = Some f }) -> f ()
  | CustomSolver _ -> 0
  | CustomSolverSens _ -> 0
  | FixedPointSolverSens _ -> c_get_num_iters rawptr
  | NewtonSolverSens _ -> c_get_num_iters rawptr
  | FixedPointSolver _ -> c_get_num_iters rawptr
  | NewtonSolver _ -> c_get_num_iters rawptr

let get_cur_iter (type d k s v) ({ rawptr; solver; _ } : (d, k, s, v) t) =
  check_compat ();
  match solver with
  | CustomSolver     (_, { get_cur_iter = Some f }) -> f ()
  | CustomSolverSens (_, { get_cur_iter = Some f }) -> f ()
  | CustomSolver _ -> 0
  | CustomSolverSens _ -> 0
  | FixedPointSolverSens _ -> c_get_cur_iter rawptr
  | NewtonSolverSens _ -> c_get_cur_iter rawptr
  | FixedPointSolver _ -> c_get_cur_iter rawptr
  | NewtonSolver _ -> c_get_cur_iter rawptr

let get_num_conv_fails (type d k s v)
                       ({ rawptr; solver; _ } : (d, k, s, v) t) =
  check_compat ();
  match solver with
  | CustomSolver     (_, { get_num_conv_fails = Some f }) -> f ()
  | CustomSolverSens (_, { get_num_conv_fails = Some f }) -> f ()
  | CustomSolver _ -> 0
  | CustomSolverSens _ -> 0
  | FixedPointSolverSens _ -> c_get_num_conv_fails rawptr
  | NewtonSolverSens _ -> c_get_num_conv_fails rawptr
  | FixedPointSolver _ -> c_get_num_conv_fails rawptr
  | NewtonSolver _ -> c_get_num_conv_fails rawptr

type ('d, 's) c_sysfn
type ('d, 's) c_lsetupfn
type ('d, 's) c_lsolvefn

type 's c_fromvaluefn

external c_call_sys_fn
  : (('d, 'k) Nvector.t, 's) c_sysfn * 's c_fromvaluefn
    -> ('d, 'k) Nvector.t
    -> ('d, 'k) Nvector.t
    -> 's
    -> unit
  = "sunml_nlsolver_call_sys_fn"

external c_call_lsetup_fn
  : (('d, 'k) Nvector.t, 's) c_lsetupfn * 's c_fromvaluefn
    -> bool
    -> 's
    -> bool
  = "sunml_nlsolver_call_lsetup_fn"

external c_call_lsolve_fn
  : (('d, 'k) Nvector.t, 's) c_lsolvefn * 's c_fromvaluefn
    -> ('d, 'k) Nvector.t
    -> 's
    -> unit
  = "sunml_nlsolver_call_lsolve_fn"

external c_call_convtest_fn' :
                 ('di, 'ki, 'si, [`Nvec]) cptr
              -> 's convtest_callback Sundials_impl.Callback.cfunptr
              -> (  ('d, 'k) Nvector.t
                  * ('d, 'k) Nvector.t
                  * float
                  * ('d, 'k) Nvector.t
                  * 's)
              -> convtest
  = "sunml_nlsolver_call_convtest_fn"

let c_call_convtest_fn cconvtestfn =
  ({ f = fun { rawptr; _ } y del tol ewt mem ->
           c_call_convtest_fn' rawptr cconvtestfn (y, del, tol, ewt, mem) }
   : 's convtest_callback)

external c_call_sys_fn_sens
  : (('d, 'k) Senswrapper.t, 's) c_sysfn * 's c_fromvaluefn
    -> ('d, 'k) Senswrapper.t
    -> ('d, 'k) Senswrapper.t
    -> 's
    -> unit
  = "sunml_nlsolver_call_sys_fn_sens"

external c_call_lsetup_fn_sens
  : (('d, 'k) Senswrapper.t, 's) c_lsetupfn * 's c_fromvaluefn
    -> bool
    -> 's
    -> bool
  = "sunml_nlsolver_call_lsetup_fn_sens"

external c_call_lsolve_fn_sens
  : (('d, 'k) Senswrapper.t, 's) c_lsolvefn * 's c_fromvaluefn
    -> ('d, 'k) Senswrapper.t
    -> 's
    -> unit
  = "sunml_nlsolver_call_lsolve_fn_sens"

external c_call_convtest_fn_sens' :
                 ('di, 'ki, 'si, [`Sens]) cptr
              -> 's convtest_callback_sens Sundials_impl.Callback.cfunptr
              -> (  ('d, 'k) Senswrapper.t
                  * ('d, 'k) Senswrapper.t
                  * float
                  * ('d, 'k) Senswrapper.t
                  * 's)
              -> convtest
  = "sunml_nlsolver_call_convtest_fn_sens"

let c_call_convtest_fn_sens cconvtestfn =
  ({ f = fun { rawptr; _ } y del tol ewt mem ->
           c_call_convtest_fn_sens' rawptr cconvtestfn (y, del, tol, ewt, mem) }
   : 's convtest_callback_sens)

module Newton = struct (* {{{ *)

  external c_make : ('d, 'k) Nvector.t
                    -> ('d, 's) callbacks Weak.t
                    -> ('d, 'k, 's, [`Nvec]) cptr
    = "sunml_nlsolver_newton_make"

  external c_make_sens
    : int
      -> ('d, 'k) Nvector.t
      -> (('d, 'k) Senswrapper.t, 's) callbacks Weak.t
      -> ('d, 'k, 's, [`Sens]) cptr
    = "sunml_nlsolver_newton_make_sens"

  (* If the returned boolean is true, then an OCaml sysfn is installed *)
  external c_get_sys_fn
    : ('d, 'k, 's, [`Nvec]) cptr
      -> bool * ((('d, 'k) Nvector.t, 's) c_sysfn * 's c_fromvaluefn) option
    = "sunml_nlsolver_newton_get_sys_fn"

  let make y =
    let callbacks = empty_callbacks () in
    {
      rawptr    = c_make y (weak_wrap callbacks);
      solver    = NewtonSolver callbacks;
      info_file = None;
      attached  = false;
    }

  let make_sens count y =
    let callbacks = empty_callbacks () in
    {
      rawptr    = c_make_sens count y (weak_wrap callbacks);
      solver    = NewtonSolverSens callbacks;
      info_file = None;
      attached  = false;
    }

  let get_sys_fn { rawptr; solver; _ } =
    check_compat ();
    (* It does not seem worth adding a phantom type argument just to avoid
       this dynamic check. *)
    match (solver : ('d, 'k, 's, [`Nvec]) solver) with
    | NewtonSolver callbacks ->
        (match c_get_sys_fn rawptr with
         | true, _ -> Some (fun y fg -> callbacks.sysfn (uw y) (uw fg))
         | _, None -> None
         | _, Some cfns -> Some (c_call_sys_fn cfns))
    | FixedPointSolver _ | CustomSolver _ -> invalid_arg "not a Newton solver"

end (* }}} *)

module FixedPoint = struct (* {{{ *)

  external c_make : ('d, 'k) Nvector.t
                    -> int
                    -> ('d, 's) callbacks Weak.t
                    -> ('d, 'k, 's, [`Nvec]) cptr
    = "sunml_nlsolver_fixedpoint_make"

  external c_make_sens
    : int
      -> ('d, 'k) Nvector.t
      -> int
      -> (('d, 'k) Senswrapper.t, 's) callbacks Weak.t
      -> ('d, 'k, 's, [`Sens]) cptr
    = "sunml_nlsolver_fixedpoint_make_sens"

  external c_get_sys_fn
    : ('d, 'k, 's, [`Nvec]) cptr
      -> bool * ((('d, 'k) Nvector.t, 's) c_sysfn * 's c_fromvaluefn) option
    = "sunml_nlsolver_fixedpoint_get_sys_fn"

  external c_set_damping : ('d, 'k, 's, 'v) cptr -> float -> unit
    = "sunml_nlsolver_fixedpoint_set_damping"

  let make ?(acceleration_vectors=0) y =
    let callbacks = empty_callbacks () in
    {
      rawptr    = c_make y acceleration_vectors (weak_wrap callbacks);
      solver    = FixedPointSolver (callbacks, acceleration_vectors);
      info_file = None;
      attached  = false;
    }

  let make_sens ?(acceleration_vectors=0) count y =
    let callbacks = empty_callbacks () in
    {
      rawptr    = c_make_sens count y acceleration_vectors (weak_wrap callbacks);
      solver    = FixedPointSolverSens (callbacks, acceleration_vectors);
      info_file = None;
      attached  = false;
    }

  let get_sys_fn { rawptr; solver; _ } =
    check_compat ();
    (* It does not seem worth adding a phantom type argument just to avoid
       this dynamic check. *)
    match (solver : ('d, 'k, 's, [`Nvec]) solver) with
    | FixedPointSolver (callbacks, _) ->
        (match c_get_sys_fn rawptr with
         | true, _ -> Some (fun y fg -> callbacks.sysfn (uw y) (uw fg))
         | _, None -> None
         | _, Some cfns -> Some (c_call_sys_fn cfns))
    | NewtonSolver _ | CustomSolver _ -> invalid_arg "not a Newton solver"

  let set_damping { rawptr; _ } beta =
    c_set_damping rawptr beta

end (* }}} *)

module Custom = struct (* {{{ *)

  external c_make
    : (('d, 'k) Nvector.t, 's) callbacks Weak.t
      -> (('d, 'k) Nvector.t, 'd, 's, [`Nvec]) ops Weak.t
      -> ('d, 'k, 's, [`Nvec]) cptr
    = "sunml_nlsolver_custom_make"

  external c_make_sens
    :     (('d, 'k) Senswrapper.t, 's) callbacks Weak.t
       -> (('d, 'k) Senswrapper.t, ('d, 'k) Senswrapper.t, 's, [`Sens]) ops Weak.t
       -> ('d, 'k, 's, [`Sens]) cptr
    = "sunml_nlsolver_custom_make_sens"

  let f_call osetf call_cfn =
    match osetf with
    | None -> ()
    | Some setf -> setf call_cfn

  (* Handle callbacks from C into a custom solver *)

  let set_c_sys_fn ops csysfn cfromvalfn =
    ops.set_sys_fn (c_call_sys_fn (csysfn, cfromvalfn))

  let set_c_lsetup_fn ops clsetupfn cfromvalfn =
    f_call ops.set_lsetup_fn (c_call_lsetup_fn (clsetupfn, cfromvalfn))

  let set_c_lsolve_fn ops clsolvefn cfromvalfn =
    f_call ops.set_lsolve_fn (c_call_lsolve_fn (clsolvefn, cfromvalfn))

  let set_c_convtest_fn ops
      (cconvtestfn : 's convtest_callback Sundials_impl.Callback.cfunptr)
    = f_call ops.set_convtest_fn
        (CConvTest Sundials_impl.Callback.{
             cptr = cconvtestfn;
             call = c_call_convtest_fn cconvtestfn })

  let set_c_sys_fn_sens ops csysfn cfromvalfn =
    ops.set_sys_fn (c_call_sys_fn_sens (csysfn, cfromvalfn))

  let set_c_lsetup_fn_sens ops clsetupfn cfromvalfn =
    f_call ops.set_lsetup_fn (c_call_lsetup_fn_sens (clsetupfn, cfromvalfn))

  let set_c_lsolve_fn_sens ops clsolvefn cfromvalfn =
    f_call ops.set_lsolve_fn (c_call_lsolve_fn_sens (clsolvefn, cfromvalfn))

  let set_c_convtest_fn_sens ops cconvtestfn =
    f_call ops.set_convtest_fn
      (CSensConvTest Sundials_impl.Callback.{
           cptr = cconvtestfn;
           call = c_call_convtest_fn_sens cconvtestfn;
         })

  let _ = Callback.register "Sundials_NonlinearSolver.set_c_sys_fn"
                            set_c_sys_fn
  let _ = Callback.register "Sundials_NonlinearSolver.set_c_lsetup_fn"
                            set_c_lsetup_fn
  let _ = Callback.register "Sundials_NonlinearSolver.set_c_lsolve_fn"
                            set_c_lsolve_fn
  let _ = Callback.register "Sundials_NonlinearSolver.set_c_convtest_fn"
                            set_c_convtest_fn

  let _ = Callback.register "Sundials_NonlinearSolver.set_c_sys_fn_sens"
                            set_c_sys_fn_sens
  let _ = Callback.register "Sundials_NonlinearSolver.set_c_lsetup_fn_sens"
                            set_c_lsetup_fn_sens
  let _ = Callback.register "Sundials_NonlinearSolver.set_c_lsolve_fn_sens"
                            set_c_lsolve_fn_sens
  let _ = Callback.register "Sundials_NonlinearSolver.set_c_convtest_fn_sens"
                            set_c_convtest_fn_sens

  (* Create custom solvers from given functions *)

  let make ?init ?setup ?set_lsetup_fn ?set_lsolve_fn ?set_convtest_fn
           ?set_max_iters ?set_info_file ?set_print_level
           ?get_num_iters ?get_cur_iter ?get_num_conv_fails
           ~nls_type ~solve ~set_sys_fn () =
    check_compat ();
    let ops = {
      nls_type;
      init;
      setup;
      solve;
      set_sys_fn;
      set_lsetup_fn;
      set_lsolve_fn;
      set_convtest_fn;
      set_max_iters;
      set_info_file;
      set_print_level;
      get_num_iters;
      get_cur_iter;
      get_num_conv_fails;
    }
    in
    let callbacks = empty_callbacks () in
    {
      rawptr    = c_make (weak_wrap callbacks) (weak_wrap ops);
      solver    = CustomSolver (callbacks, ops);
      info_file = None;
      attached  = false;
    }

  let make_sens
                ?init ?setup ?set_lsetup_fn ?set_lsolve_fn ?set_convtest_fn
                ?set_max_iters ?set_info_file ?set_print_level
                ?get_num_iters ?get_cur_iter ?get_num_conv_fails
                ~nls_type ~solve ~set_sys_fn () =
    check_compat ();
    let ops = {
      nls_type;
      init;
      setup;
      solve;
      set_sys_fn;
      set_lsetup_fn;
      set_lsolve_fn;
      set_convtest_fn;
      set_max_iters;
      set_info_file;
      set_print_level;
      get_num_iters;
      get_cur_iter;
      get_num_conv_fails;
    }
    in
    let callbacks = empty_callbacks () in
    {
      rawptr    = c_make_sens (weak_wrap callbacks) (weak_wrap ops);
      solver    = CustomSolverSens (callbacks, ops);
      info_file = None;
      attached  = false;
    }

end (* }}} *)

(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "sunml_nlsolver_init_module"

let _ =
  c_init_module
    (* Exceptions must be listed in the same order as lsolver_exn_index.  *)
    [|
      VectorOpError;
      IncorrectUse;
      ExtFail;
    |]

