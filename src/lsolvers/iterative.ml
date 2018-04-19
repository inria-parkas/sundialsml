(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2018 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(* "Simulate" Linear Solvers in Sundials < 3.0.0 *)
let in_compat_mode =
  match Sundials.sundials_version with
  | 2,_,_ -> true
  | _ -> false

include Lsolver_impl.Iterative

external c_set_maxl
  : ('nd, 'nk) cptr
    -> ('nd, 'nk, [< `Spbcgs|`Sptfqmr|`Pcg]) solver
    -> int
    -> unit
  = "ml_lsolver_set_maxl"

external c_set_gs_type
  : ('nd, 'nk) cptr
    -> ('nd, 'nk, [< `Spfgmr|`Spgmr]) solver
    -> gramschmidt_type
    -> unit
  = "ml_lsolver_set_gs_type"

external c_set_max_restarts
  : ('nd, 'nk) cptr
    -> ('nd, 'nk, [< `Spfgmr|`Spgmr]) solver
    -> int
    -> unit
  = "ml_lsolver_set_max_restarts"

let set_maxl { rawptr; solver; compat } maxl =
  if in_compat_mode then compat.set_maxl maxl
  else c_set_maxl rawptr solver maxl

let set_gs_type { rawptr; solver; compat } gs_type =
  if in_compat_mode then compat.set_gs_type gs_type
  else c_set_gs_type rawptr solver gs_type

let set_max_restarts { rawptr; solver; compat } max_restarts =
  if in_compat_mode then raise Sundials.NotImplementedBySundialsVersion
  else c_set_max_restarts rawptr solver max_restarts

let set_prec_type { rawptr; solver; compat; check_prec_type } prec_type =
  if not (check_prec_type prec_type) then raise Lsolver.IllegalPrecType;
  if in_compat_mode then compat.set_prec_type prec_type
  else c_set_prec_type rawptr solver prec_type

let default = function
  | Some x -> x
  | None -> 0

external c_spbcgs : int -> ('d, 'k) Nvector.t -> ('nd, 'nk) cptr
  = "ml_lsolver_spbcgs"

let spbcgs ?maxl nvec =
  let maxl = default maxl in
  let cptr = c_spbcgs maxl nvec in
  let compat =
    if in_compat_mode
    then let r = { info with maxl = maxl } in
         r.set_maxl <- (fun m -> r.maxl <- m);
         r
    else info
  in
  { rawptr = cptr;
    solver = Spbcgs;
    compat = compat;
    check_prec_type = (fun _ -> true);
    attached = false;
  }

external c_spfgmr : int -> ('d, 'k) Nvector.t -> ('nd, 'nk) cptr
  = "ml_lsolver_spfgmr"

let spfgmr ?maxl ?max_restarts ?gs_type nvec =
  let maxl = default maxl in
  let cptr = c_spfgmr maxl nvec in
  let compat =
    if in_compat_mode
    then begin
           if max_restarts <> None
             then raise Sundials.NotImplementedBySundialsVersion;
           let r = { info with maxl = maxl; gs_type = gs_type } in
           r.set_gs_type <- (fun t -> r.gs_type <- Some t); r
         end
    else begin
      (match max_restarts with
       | Some mr -> c_set_max_restarts cptr Spfgmr mr
       | _ -> ());
      (match gs_type with
       | Some gst -> c_set_gs_type cptr Spfgmr gst
       | _ -> ());
      info
    end
  in
  { rawptr = cptr;
    solver = Spfgmr;
    compat = compat;
    check_prec_type = (fun _ -> true);
    attached = false;
  }

external c_spgmr : int -> ('d, 'k) Nvector.t -> ('nd, 'nk) cptr
  = "ml_lsolver_spgmr"

let spgmr ?maxl ?max_restarts ?gs_type nvec =
  let maxl = default maxl in
  let cptr = c_spgmr maxl nvec in
  let compat =
    if in_compat_mode
    then begin
           if max_restarts <> None
             then raise Sundials.NotImplementedBySundialsVersion;
           let r = { info with maxl = maxl; gs_type = gs_type } in
           r.set_gs_type <- (fun t -> r.gs_type <- Some t);
           r
         end
    else begin
      (match max_restarts with
       | Some mr -> c_set_max_restarts cptr Spgmr mr
       | _ -> ());
      (match gs_type with
       | Some gst -> c_set_gs_type cptr Spgmr gst
       | _ -> ());
      info
    end
  in
  { rawptr = cptr;
    solver = Spgmr;
    compat = compat;
    check_prec_type = (fun _ -> true);
    attached = false;
  }

external c_sptfqmr : int -> ('d, 'k) Nvector.t -> ('nd, 'nk) cptr
  = "ml_lsolver_sptfqmr"

let sptfqmr ?maxl nvec =
  let maxl = default maxl in
  let cptr = c_sptfqmr maxl nvec in
  let compat =
    if in_compat_mode
    then let r = { info with maxl = maxl } in
         r.set_maxl <- (fun m -> r.maxl <- m);
         r
    else info
  in
  { rawptr = cptr;
    solver = Sptfqmr;
    compat = compat;
    check_prec_type = (fun _ -> true);
    attached = false;
  }

external c_pcg : int -> ('d, 'k) Nvector.t -> ('nd, 'nk) cptr
  = "ml_lsolver_pcg"

let pcg ?maxl nvec =
  let maxl = default maxl in
  let cptr = c_pcg maxl nvec in
  let compat =
    if in_compat_mode
    then let r = { info with maxl = maxl } in
         r.set_maxl <- (fun m -> r.maxl <- m); r
    else info
  in
  { rawptr = cptr;
    solver = Pcg;
    compat = compat;
    check_prec_type = (fun _ -> true);
    attached = false;
  }

module Custom = struct (* {{{ *)

  type 'lsolver tag = [`Custom of 'lsolver]

  type ('data, 'kind) atimesfn =
       ('data, 'kind) Nvector.t
    -> ('data, 'kind) Nvector.t
    -> unit

  type psetupfn = unit -> unit

  type ('data, 'kind) psolvefn =
       ('data, 'kind) Nvector.t
    -> ('data, 'kind) Nvector.t
    -> float
    -> bool
    -> unit

  type ('data, 'kind, 'lsolver) ops = {

    init : 'lsolver -> unit;

    setup : 'lsolver -> unit;

    solve : 'lsolver
            -> ('data, 'kind) Nvector.t
            -> ('data, 'kind) Nvector.t
            -> float
            -> unit;

    set_atimes
      : ('lsolver -> ('data, 'kind) atimesfn -> unit) option;

    set_preconditioner
      : ('lsolver 
         -> psetupfn option
         -> ('data, 'kind) psolvefn option
         -> unit) option;

    set_scaling_vectors
      : ('lsolver
         -> ('data, 'kind) Nvector.t option
         -> ('data, 'kind) Nvector.t option
         -> unit) option;

    get_num_iters : ('lsolver -> int) option;

    get_res_norm : ('lsolver -> float) option;

    get_res_id : ('lsolver -> ('data, 'kind) Nvector.t) option;

    get_work_space : ('lsolver -> int * int) option;
  }

  let wrap_set_atimes fseto ldata =
    match fseto with
    | None ->
        fun _  -> failwith ("internal error: Iterative.Custom.set_atimes")
    | Some fset ->
        let fset' = fset ldata in
        fun fd -> fset' (Lsolver_impl.Custom.call_atimes fd)

  let wrap_set_preconditioner fseto ldata =
    match fseto with
    | None -> fun _  ->
        failwith ("internal error: Iterative.Custom.set_preconditioner")
    | Some fset ->
        let fset' = fset ldata in
        fun fd has_setup has_solve -> fset'
          (if has_setup
           then Some (fun () -> Lsolver_impl.Custom.call_psetup fd) else None)
          (if has_solve
           then Some (Lsolver_impl.Custom.call_psolve fd) else None)

  let mapo s fo x =
    match fo with
    | None -> (fun _ -> failwith ("internal error: Iterative.Custom." ^ s))
    | Some f -> f x

  let mapu s fo x =
    match fo with
    | None -> (fun _ -> failwith ("internal error: Iterative.Custom." ^ s))
    | Some f -> fun () -> f x

  let make { init = finit;
             setup = fsetup;
             solve = fsolve;
             set_atimes = fset_atimes;
             set_preconditioner = fset_preconditioner;
             set_scaling_vectors = fset_scaling_vectors;
             get_num_iters = fget_num_iters;
             get_res_norm = fget_res_norm;
             get_res_id = fget_res_id;
             get_work_space = fget_work_space } ldata =
    match Sundials.sundials_version with
    | 2,_,_ -> raise Sundials.NotImplementedBySundialsVersion;
    | _ -> ();
    let ops = Lsolver_impl.Custom.({
            init = (fun () -> finit ldata);
            setup = (fun () -> fsetup ldata);
            solve = (fun () -> fsolve ldata);
            set_atimes = wrap_set_atimes fset_atimes ldata;
            set_preconditioner =
              wrap_set_preconditioner fset_preconditioner ldata;
            set_scaling_vectors =
              mapo "set_scaling_vectors" fset_scaling_vectors ldata;
            get_num_iters = mapu "get_num_iters" fget_num_iters ldata;
            get_res_norm = mapu "get_res_norm" fget_res_norm ldata;
            get_res_id = mapu "get_res_id" fget_res_id ldata;
            get_work_space = mapu "get_work_space" fget_work_space ldata;
          }) in
    let only_ops = Lsolver_impl.Custom.({
      has_set_atimes          = fset_atimes <> None;
      has_set_preconditioner  = fset_preconditioner <> None;
      has_set_scaling_vectors = fset_scaling_vectors <> None;
      has_get_num_iters       = fget_num_iters <> None;
      has_get_res_norm        = fget_res_norm <> None;
      has_get_res_id          = fget_res_id <> None;
      has_get_work_space      = fget_work_space <> None;
    })
    in Lsolver_impl.Iterative.({
        rawptr = c_make_custom 1 ops only_ops;
        solver = Custom (ldata, ops);
        compat = info;
        check_prec_type = (fun _ -> true);
        attached = false;
      })

  let unwrap { Lsolver_impl.Iterative.solver = Custom (ldata, _) } = ldata

end (* }}} *)

module Algorithms = struct (* {{{ *)

  external qr_fact : Sundials.RealArray2.t
                     -> Sundials.RealArray.t
                     -> bool
                     -> unit
      = "c_spils_qr_fact"
   
  external qr_sol : Sundials.RealArray2.t
                   -> Sundials.RealArray.t
                   -> Sundials.RealArray.t
                   -> unit
      = "c_spils_qr_sol"

  external modified_gs : (('a, 'k) Nvector.t) array
                         -> Sundials.RealArray2.t
                         -> int
                         -> int
                         -> float
      = "c_spils_modified_gs"

  external classical_gs' : (('a, 'k) Nvector.t) array
                           * Sundials.RealArray2.t
                           * int
                           * int
                           * ('a, 'k) Nvector.t
                           * Sundials.RealArray.t
                           -> float
      = "c_spils_classical_gs"

  let classical_gs v h k p temp s = classical_gs' (v, h, k, p, temp, s)

end (* }}} *)

