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
open Sundials_impl

include Sundials_NonlinearSolver_impl

type ('data, 'kind, 's) t
    = ('data, 'kind, 's) Sundials_NonlinearSolver_impl.nonlinear_solver

let check_compat () =
  if Versions.in_compat_mode2_3
  then raise Config.NotImplementedBySundialsVersion

external c_init : ('d, 'k, 's) cptr -> unit
  = "sunml_nlsolver_init"

external c_setup : ('d, 'k, 's) cptr -> ('d, 'k) Nvector.t -> unit
  = "sunml_nlsolver_setup"

external c_solve :
  ('d, 'k, 's) cptr
  -> (('d, 'k) Nvector.t * ('d, 'k) Nvector.t * ('d, 'k) Nvector.t
      * float * bool)
  -> unit
  = "sunml_nlsolver_solve"

external c_set_sys_fn      : ('d, 'k, 's) cptr -> unit
  = "sunml_nlsolver_set_sys_fn"

external c_set_lsetup_fn   : ('d, 'k, 's) cptr -> unit
  = "sunml_nlsolver_set_lsetup_fn"

external c_set_lsolve_fn   : ('d, 'k, 's) cptr -> unit
  = "sunml_nlsolver_set_lsolve_fn"

external c_set_convtest_fn : ('d, 'k, 's) cptr -> unit
  = "sunml_nlsolver_set_convtest_fn"

external c_set_max_iters : ('d, 'k, 's) cptr -> int -> unit
  = "sunml_nlsolver_set_max_iters"

external c_get_num_iters : ('d, 'k, 's) cptr -> int
  = "sunml_nlsolver_get_num_iters"

external c_get_cur_iter  : ('d, 'k, 's) cptr -> int
  = "sunml_nlsolver_get_cur_iter"

external c_get_num_conv_fails : ('d, 'k, 's) cptr -> int
  = "sunml_nlsolver_get_num_conv_fails"

external c_set_info_file_fixedpoint : ('d, 'k, 's) cptr -> Logfile.t -> unit
  = "sunml_nlsolver_set_info_file_fixedpoint"

external c_set_info_file_newton : ('d, 'k, 's) cptr -> Logfile.t -> unit
  = "sunml_nlsolver_set_info_file_newton"

external c_set_print_level_fixedpoint : ('d, 'k, 's) cptr -> int -> unit
  = "sunml_nlsolver_set_print_level_fixedpoint"

external c_set_print_level_newton : ('d, 'k, 's) cptr -> int -> unit
  = "sunml_nlsolver_set_print_level_newton"

(* - - - OCaml invoking init/setup/solve - - - *)

let init (type d k s) ({ rawptr; solver } : (d, k, s) t) =
  check_compat ();
  match solver with
  | CustomSolver     { init = Some f }  -> f ()              (* O/Onls *)
  | CustomSensSolver { init = Some f }  -> f ()
  | CustomSolver _     -> ()
  | CustomSensSolver _ -> ()
  | FixedPointSolver _ | NewtonSolver     -> c_init rawptr   (* O/Cnls *)

let setup { rawptr; solver } ~y =
  check_compat ();
  match solver with
  | CustomSolver     { setup = Some f } -> f y ()            (* O/Onls *)
  | CustomSolver     { setup = None   } -> ()
  | FixedPointSolver _ | NewtonSolver
      -> c_setup rawptr y                                    (* O/Cnls *)

let solve { rawptr; solver } ~y0 ~ycor ~w tol callLSetup =
  check_compat ();
  match solver with
  | CustomSolver { solve = f } -> f y0 ycor w tol callLSetup () (* O/Onls *)
  | FixedPointSolver _ | NewtonSolver
      -> c_solve rawptr (y0, ycor, w, tol, callLSetup)          (* O/Cnls *)

(* - - - OCaml callback configuration - - - *)

let set_sys_fn { rawptr; solver; callbacks } cbf =
  check_compat ();
  match (solver : ('d, 'k, user) solver) with
  | CustomSolver { set_sys_fn = set }                        (* O/Onls *)
      -> set (fun y fg -> Nvector.(cbf (unwrap y) (unwrap fg)))
  | FixedPointSolver _ | NewtonSolver                        (* O/Cnls *)
      -> callbacks.sysfn <- cbf;
         c_set_sys_fn rawptr

let set_lsetup_fn { rawptr; solver; callbacks } cbf =
  check_compat ();
  match (solver : ('d, 'k, user) solver) with
  | CustomSolver { set_lsetup_fn = Some set } -> set cbf     (* O/Onls *)
  | CustomSolver { set_lsetup_fn = None } -> ()
  | FixedPointSolver _ | NewtonSolver                        (* O/Cnls *)
      -> callbacks.lsetupfn <- cbf;
         c_set_lsetup_fn rawptr

let set_lsolve_fn { rawptr; solver; callbacks } cbf =
  check_compat ();
  match (solver : ('d, 'k, user) solver) with
  | CustomSolver { set_lsolve_fn = Some set }                (* O/Onls *)
      -> set (fun b -> Nvector.(cbf (unwrap b)))
  | CustomSolver { set_lsolve_fn = None } -> ()
  | FixedPointSolver _ | NewtonSolver                        (* O/Cnls *)
      -> callbacks.lsolvefn <- cbf;
         c_set_lsolve_fn rawptr

let set_convtest_fn (type d k s)
                    ({ rawptr; solver; callbacks } : (d, k, s) t)
                    (cbf : (d, s) convtestfn) =
  check_compat ();
  match solver with
  | CustomSolver     { set_convtest_fn = Some set } ->       (* O/Onls *)
      set (fun y del tol ewt ->
             Nvector.(cbf (unwrap y) (unwrap del) tol (unwrap ewt)))
  | CustomSensSolver { set_convtest_fn = Some set } -> set cbf
  | CustomSolver _ -> ()
  | CustomSensSolver _ -> ()
  | FixedPointSolver _ | NewtonSolver                        (* O/Cnls *)
      -> callbacks.convtestfn <- cbf;
         c_set_convtest_fn rawptr

let set_max_iters (type d k s) ({ rawptr; solver } : (d, k, s) t) i =
  check_compat ();
  match solver with
  | CustomSolver     { set_max_iters = Some f } -> f i
  | CustomSensSolver { set_max_iters = Some f } -> f i
  | CustomSolver _ -> ()
  | CustomSensSolver _ -> ()
  | FixedPointSolver _ | NewtonSolver
      -> c_set_max_iters rawptr i

let set_info_file (type d k s) ({ rawptr; solver } : (d, k, s) t) file =
  if Versions.sundials_lt530 then raise Config.NotImplementedBySundialsVersion;
  match solver with
  | CustomSolver     { set_info_file = Some f } -> f file
  | CustomSensSolver { set_info_file = Some f } -> f file
  | CustomSolver _ -> ()
  | CustomSensSolver _ -> ()
  | FixedPointSolver _ -> c_set_info_file_fixedpoint rawptr file
  | NewtonSolver -> c_set_info_file_newton rawptr file

let set_print_level (type d k s) ({ rawptr; solver } : (d, k, s) t) level =
  if Versions.sundials_lt530 then raise Config.NotImplementedBySundialsVersion;
  let level = if level then 1 else 0 in
  match solver with
  | CustomSolver     { set_print_level = Some f } -> f level
  | CustomSensSolver { set_print_level = Some f } -> f level
  | CustomSolver _ -> ()
  | CustomSensSolver _ -> ()
  | FixedPointSolver _ -> c_set_print_level_fixedpoint rawptr level
  | NewtonSolver -> c_set_print_level_newton rawptr level

let get_num_iters (type d k s) ({ rawptr; solver } : (d, k, s) t) =
  check_compat ();
  match solver with
  | CustomSolver     { get_num_iters = Some f } -> f ()
  | CustomSensSolver { get_num_iters = Some f } -> f ()
  | CustomSolver _     -> 0
  | CustomSensSolver _ -> 0
  | FixedPointSolver _ | NewtonSolver
      -> c_get_num_iters rawptr

let get_cur_iter (type d k s) ({ rawptr; solver } : (d, k, s) t) =
  check_compat ();
  match solver with
  | CustomSolver     { get_cur_iter = Some f } -> f ()
  | CustomSensSolver { get_cur_iter = Some f } -> f ()
  | CustomSolver _     -> 0
  | CustomSensSolver _ -> 0
  | FixedPointSolver _ | NewtonSolver
      -> c_get_cur_iter rawptr

let get_num_conv_fails (type d k s)
                       ({ rawptr; solver } : (d, k, s) t) =
  check_compat ();
  match solver with
  | CustomSolver     { get_num_conv_fails = Some f } -> f ()
  | CustomSensSolver { get_num_conv_fails = Some f } -> f ()
  | CustomSolver _     -> 0
  | CustomSensSolver _ -> 0
  | FixedPointSolver _ | NewtonSolver
      -> c_get_num_conv_fails rawptr

type ('d, 's) c_sysfn
type ('d, 's) c_lsetupfn
type ('d, 's) c_lsolvefn
type ('d, 's) c_convtestfn

external c_call_sys_fn : (('d, 'k) Nvector.t, 's) c_sysfn
              -> ('d, 'k) Nvector.t -> ('d, 'k) Nvector.t -> 's -> unit
  = "sunml_nlsolver_call_sys_fn"

external c_call_lsetup_fn : (('d, 'k) Nvector.t, 's) c_lsetupfn
              -> bool -> 's -> bool
  = "sunml_nlsolver_call_lsetup_fn"

external c_call_lsolve_fn : (('d, 'k) Nvector.t, 's) c_lsolvefn
              -> ('d, 'k) Nvector.t -> 's -> unit
  = "sunml_nlsolver_call_lsolve_fn"

external c_call_convtest_fn' :
              ('d, 'k, 's) cptr
              -> (('d, 'k) Nvector.t, 's) c_convtestfn
              -> (  ('d, 'k) Nvector.t
                  * ('d, 'k) Nvector.t
                  * float
                  * ('d, 'k) Nvector.t
                  * 's)
              -> convtest
  = "sunml_nlsolver_call_convtest_fn"

let c_call_convtest_fn nlsptr cfn y del tol ewt mem =
  c_call_convtest_fn' nlsptr cfn (y, del, tol, ewt, mem)

external c_call_sys_fn_sens : (('d, 'k) Senswrapper.t, 's) c_sysfn
              -> ('d, 'k) Senswrapper.t -> ('d, 'k) Senswrapper.t -> 's -> unit
  = "sunml_nlsolver_call_sys_fn_sens"

external c_call_lsetup_fn_sens : (('d, 'k) Senswrapper.t, 's) c_lsetupfn
              -> bool -> 's -> bool
  = "sunml_nlsolver_call_lsetup_fn_sens"

external c_call_lsolve_fn_sens : (('d, 'k) Senswrapper.t, 's) c_lsolvefn
              -> ('d, 'k) Senswrapper.t -> 's -> unit
  = "sunml_nlsolver_call_lsolve_fn_sens"

external c_call_convtest_fn_sens' :
              ('d, 'k, 's) cptr
              -> (('d, 'k) Senswrapper.t, 's) c_convtestfn
              -> (  ('d, 'k) Senswrapper.t
                  * ('d, 'k) Senswrapper.t
                  * float
                  * ('d, 'k) Senswrapper.t
                  * 's)
              -> convtest
  = "sunml_nlsolver_call_convtest_fn_sens"

let c_call_convtest_fn_sens nlsptr cfn y del tol ewt mem =
  c_call_convtest_fn_sens' nlsptr cfn (y, del, tol, ewt, mem)

module Newton = struct (* {{{ *)

  external c_make
    : ('d, 'k) Nvector.t -> ('d, 'k, 's) callbacks -> ('d, 'k, 's) cptr
    = "sunml_nlsolver_newton_make"

  external c_make_sens
    : int
      -> ('d, 'k) Nvector.t
      -> (('d, 'k) Senswrapper.t, 'k, 'a integrator) callbacks
      -> (('d, 'k) Senswrapper.t, 'k, 'a integrator) cptr
    = "sunml_nlsolver_newton_make_sens"

  external c_get_sys_fn
    : ('d, 'k, 's) cptr -> (('d, 'k) Nvector.t, 's) c_sysfn option
    = "sunml_nlsolver_newton_get_sys_fn"

  let make y =
    let callbacks = empty_callbacks () in
    {
      rawptr    = c_make y callbacks;
      solver    = NewtonSolver;
      callbacks = callbacks;
      attached  = false;
    }

  let make_sens count y =
    let callbacks = empty_callbacks () in
    {
      rawptr    = c_make_sens count y callbacks;
      solver    = NewtonSolver;
      callbacks = callbacks;
      attached  = false;
    }

  let get_sys_fn { rawptr; solver } =
    check_compat ();
    (* It does not seem worth adding a phantom type argument just to avoid
       this dynamic check. *)
    if solver <> NewtonSolver then invalid_arg "not a Newton solver";
    match c_get_sys_fn rawptr with
    | None -> None
    | Some f -> Some (c_call_sys_fn f)

end (* }}} *)

module FixedPoint = struct (* {{{ *)

  external c_make
    : ('d, 'k) Nvector.t -> int -> ('d, 'k, 's) callbacks -> ('d, 'k, 's) cptr
    = "sunml_nlsolver_fixedpoint_make"

  external c_make_sens
    : int
      -> ('d, 'k) Nvector.t
      -> int
      -> (('d, 'k) Senswrapper.t, 'k, 'a integrator) callbacks
      -> (('d, 'k) Senswrapper.t, 'k, 'a integrator) cptr
    = "sunml_nlsolver_fixedpoint_make_sens"

  external c_get_sys_fn
    : ('d, 'k, 's) cptr -> (('d, 'k) Nvector.t, 's) c_sysfn option
    = "sunml_nlsolver_fixedpoint_get_sys_fn"

  external c_set_damping : ('d, 'k, 's) cptr -> float -> unit
    = "sunml_nlsolver_fixedpoint_set_damping"

  let make ?(acceleration_vectors=0) y =
    let callbacks = empty_callbacks () in
    {
      rawptr    = c_make y acceleration_vectors callbacks;
      solver    = FixedPointSolver acceleration_vectors;
      callbacks = callbacks;
      attached  = false;
    }

  let make_sens ?(acceleration_vectors=0) count y =
    let callbacks = empty_callbacks () in
    {
      rawptr    = c_make_sens count y acceleration_vectors callbacks;
      solver    = FixedPointSolver acceleration_vectors;
      callbacks = callbacks;
      attached  = false;
    }

  (* It does not seem worth adding a phantom type argument just to avoid
     this dynamic check. *)
  let get_sys_fn { rawptr; solver } =
    check_compat ();
    (match solver with
     | FixedPointSolver _ -> ()
     | _ -> invalid_arg "not a FixedPoint solver");
    match c_get_sys_fn rawptr with
    | None -> None
    | Some f -> Some (c_call_sys_fn f)

  let set_damping { rawptr; _ } beta =
    c_set_damping rawptr beta

end (* }}} *)

module Custom = struct (* {{{ *)

  external c_make
    : ('d, 'k, 's) callbacks -> ('nd, 's) ops -> ('d, 'k, 's) cptr
    = "sunml_nlsolver_custom_make"

  external c_make_sens
    :     (('d, 'k) Senswrapper.t, 'k, 's) callbacks
       -> (('d, 'k) Senswrapper.t, 's) ops
       -> (('d, 'k) Senswrapper.t, 'k, 's) cptr
    = "sunml_nlsolver_custom_make_sens"

  let f_call osetf call_cfn cfn =
    match osetf with
    | None -> ()
    | Some setf -> setf (call_cfn cfn)

  (* Handle callbacks from C into a custom solver *)

  let set_c_sys_fn ops csysfn = ops.set_sys_fn (c_call_sys_fn csysfn)
  let set_c_lsetup_fn ops     = f_call ops.set_lsetup_fn c_call_lsetup_fn
  let set_c_lsolve_fn ops     = f_call ops.set_lsolve_fn c_call_lsolve_fn
  let set_c_convtest_fn nlsptr ops = f_call ops.set_convtest_fn
                                            (c_call_convtest_fn nlsptr)

  let set_c_sys_fn_sens ops csysfn = ops.set_sys_fn (c_call_sys_fn_sens csysfn)
  let set_c_lsetup_fn_sens ops     = f_call ops.set_lsetup_fn c_call_lsetup_fn_sens
  let set_c_lsolve_fn_sens ops     = f_call ops.set_lsolve_fn c_call_lsolve_fn_sens
  let set_c_convtest_fn_sens nlsptr ops = f_call ops.set_convtest_fn
                                            (c_call_convtest_fn_sens nlsptr)

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
      rawptr    = c_make callbacks ops;
      solver    = CustomSolver ops;

      (* These callbacks are never modified or used since calls to sysfn,
         lsetupfn, lsolvefn, and convtestfn occur directly from OCaml. *)
      callbacks = callbacks;
      attached  = false;
    }

  let make_sens
                ?init ?setup ?set_lsetup_fn ?set_lsolve_fn ?set_convtest_fn
                ?set_max_iters ?set_info_file ?set_print_level
                ?get_num_iters ?get_cur_iter ?get_num_conv_fails
                ~nls_type ~solve ~set_sys_fn () =
    check_compat ();
    let ops = ({
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
    } : (('d, 'k) Senswrapper.t, 'a integrator) ops)
    in
    let callbacks = empty_callbacks () in
    {
      rawptr    = c_make_sens callbacks ops;
      solver    = CustomSensSolver ops;

      (* These callbacks are never modified or used since calls to sysfn,
         lsetupfn, lsolvefn, and convtestfn occur directly from OCaml. *)
      callbacks = callbacks;
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

