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

open Sundials
module LSI = Sundials_LinearSolver_impl

include LSI

exception InvalidLinearSolver
exception UnrecoverableFailure of bool
exception MatrixNotSquare
exception MatrixVectorMismatch
exception InsufficientStorageUpperBandwidth
exception ATimesFailure of bool
exception PSetFailure of bool
exception PSolveFailure of bool
exception GSFailure
exception QRSolFailure
exception VectorOpError
exception ResReduced
exception ConvFailure
exception QRfactFailure
exception LUfactFailure
exception PackageFailure of bool
exception IllegalPrecType
exception InternalFailure of (string * int)
exception ZeroInDiagonal of int

(* "Simulate" Linear Solvers in Sundials < 3.0.0 *)
let in_compat_mode =
  match Config.sundials_version with
  | 2,_,_ -> true
  | _ -> false

type ('mat, 'data, 'kind, 'tag) t = ('mat, 'data, 'kind, 'tag) linear_solver

(** Alias for linear solvers that are restricted to serial nvectors. *)
type ('mat, 'kind, 'tag) serial_t
  = ('mat, Nvector_serial.data, [>Nvector_serial.kind] as 'kind, 'tag)
      linear_solver

module Direct = struct (* {{{ *)

  external c_dense
           : 'k Nvector_serial.any
             -> 'k Matrix.dense
             -> (Matrix.Dense.t, Nvector_serial.data, 'k) cptr
    = "sunml_lsolver_dense"

  let dense nvec mat = LS {
    rawptr = c_dense nvec mat;
    solver = Dense;
    matrix = Some mat;
    compat = LSI.Iterative.info;
    check_prec_type = (fun _ -> true);
    attached = false;
  }

  external c_lapack_dense
           : 'k Nvector_serial.any
             -> 'k Matrix.dense
             -> (Matrix.Dense.t, Nvector_serial.data, 'k) cptr
    = "sunml_lsolver_lapack_dense"

  let lapack_dense nvec mat =
    if not Config.lapack_enabled
    then raise Config.NotImplementedBySundialsVersion;
    LS {
      rawptr = c_lapack_dense nvec mat;
      solver = LapackDense;
      matrix = Some mat;
      compat = LSI.Iterative.info;
      check_prec_type = (fun _ -> true);
      attached = false;
    }

  external c_band
           : 'k Nvector_serial.any
             -> 'k Matrix.band
             -> (Matrix.Band.t, Nvector_serial.data, 'k) cptr
    = "sunml_lsolver_band"

  let band nvec mat = LS {
    rawptr = c_band nvec mat;
    solver = Band;
    matrix = Some mat;
    compat = LSI.Iterative.info;
    check_prec_type = (fun _ -> true);
    attached = false;
  }

  external c_lapack_band
           : 'k Nvector_serial.any
             -> 'k Matrix.band
             -> (Matrix.Band.t, Nvector_serial.data, 'k) cptr
    = "sunml_lsolver_lapack_band"

  let lapack_band nvec mat =
    if not Config.lapack_enabled
    then raise Config.NotImplementedBySundialsVersion;
    LS {
      rawptr = c_lapack_band nvec mat;
      solver = LapackBand;
      matrix = Some mat;
      compat = LSI.Iterative.info;
      check_prec_type = (fun _ -> true);
      attached = false;
    }

  module Klu = struct (* {{{ *)
    include LSI.Klu

    external c_klu
             : 'k Nvector_serial.any
               -> ('s, 'k) Matrix.sparse
               -> ('s Matrix.Sparse.t, Nvector_serial.data, 'k) cptr
      = "sunml_lsolver_klu"

    let make ?ordering nvec mat =
      if not Config.klu_enabled
      then raise Config.NotImplementedBySundialsVersion;
      let cptr = c_klu nvec mat in
      let info =
        if in_compat_mode
        then let r = { (info()) with ordering = ordering } in
             r.set_ordering <- (fun o -> r.ordering <- Some o); r
        else info ()
      in
      LS {
        rawptr = cptr;
        solver = Klu info;
        matrix = Some mat;
        compat = LSI.Iterative.info;
        check_prec_type = (fun _ -> true);
        attached = false;
      }

    external c_reinit
             : ('s Matrix.Sparse.t, Nvector_serial.data, 'k) cptr
               -> ('s, 'k) Matrix.sparse
               -> unit
      = "sunml_lsolver_klu_reinit"

    let reinit (LS { rawptr = cptr; solver }) mat ?nnz () =
      if in_compat_mode then
        match solver with
        | Klu { reinit = f } ->
           let smat = Matrix.unwrap mat in
           f (fst (Matrix.Sparse.size smat)) nnz
        | _ -> assert false
      else begin
          (match nnz with
           | Some n -> Matrix.Sparse.resize ~nnz:n (Matrix.unwrap mat)
           | None -> ());
          c_reinit cptr mat
        end

    external c_set_ordering
             : ('s Matrix.Sparse.t, Nvector_serial.data, 'k) cptr
               -> ordering
               -> unit
      = "sunml_lsolver_klu_set_ordering"

    let set_ordering (LS { rawptr = cptr; solver }) ordering =
      if in_compat_mode then
        match solver with
        | Klu { set_ordering = f } -> f ordering
        | _ -> assert false
      else c_set_ordering cptr ordering

  end (* }}} *)

  let klu = Klu.make

  module Superlumt = struct (* {{{ *)
    include LSI.Superlumt

    external c_superlumt
             : 'k Nvector_serial.any
               -> ('s, 'k) Matrix.sparse
               -> int
               -> ('s Matrix.Sparse.t, Nvector_serial.data, 'k) cptr
      = "sunml_lsolver_superlumt"

    let make ?ordering ~nthreads nvec mat =
      if not Config.superlumt_enabled
         || (in_compat_mode && not Matrix.(Sparse.is_csc (unwrap mat)))
      then raise Config.NotImplementedBySundialsVersion;
      let cptr = c_superlumt nvec mat nthreads in
      let info =
        if in_compat_mode
        then let r = { (info nthreads) with ordering = ordering } in
             r.set_ordering <- (fun o -> r.ordering <- Some o); r
        else info nthreads
      in
      LS {
        rawptr = cptr;
        solver = Superlumt info;
        matrix = Some mat;
        compat = LSI.Iterative.info;
        check_prec_type = (fun _ -> true);
        attached = false;
      }

    external c_set_ordering
             : ('s Matrix.Sparse.t, Nvector_serial.data, 'k) cptr
               -> ordering
               -> unit
      = "sunml_lsolver_superlumt_set_ordering"

    let set_ordering (LS { rawptr = cptr; solver }) ordering =
      if in_compat_mode then
        match solver with
        | Superlumt { set_ordering = f } -> f ordering
        | _ -> assert false
      else c_set_ordering cptr ordering

  end (* }}} *)

  let superlumt = Superlumt.make

end (* }}} *)

module Iterative = struct (* {{{ *)
  include LSI.Iterative

  external c_set_maxl
           : ('m, 'nd, 'nk) cptr
             -> ('m, 'nd, 'nk, [< `Iter|`Spbcgs|`Sptfqmr|`Pcg]) solver_data
             -> int
             -> unit
    = "sunml_lsolver_set_maxl"

  external c_set_gs_type
           : ('m, 'nd, 'nk) cptr
             -> ('m, 'nd, 'nk, [< `Iter|`Spfgmr|`Spgmr]) solver_data
             -> gramschmidt_type
             -> unit
    = "sunml_lsolver_set_gs_type"

  external c_set_max_restarts
           : ('m, 'nd, 'nk) cptr
             -> ('m, 'nd, 'nk, [< `Iter|`Spfgmr|`Spgmr]) solver_data
             -> int
             -> unit
    = "sunml_lsolver_set_max_restarts"

  let set_maxl (LS { rawptr; solver; compat }) maxl =
    if in_compat_mode then compat.set_maxl maxl
    else c_set_maxl rawptr solver maxl

  let set_gs_type (LS { rawptr; solver; compat }) gs_type =
    if in_compat_mode then compat.set_gs_type gs_type
    else c_set_gs_type rawptr solver gs_type

  let set_max_restarts (LS { rawptr; solver; compat }) max_restarts =
    if in_compat_mode then compat.set_max_restarts max_restarts
    else c_set_max_restarts rawptr solver max_restarts

  let set_prec_type (LS { rawptr; solver; compat; check_prec_type }) prec_type =
    if not (check_prec_type prec_type) then raise IllegalPrecType;
    if in_compat_mode then compat.set_prec_type prec_type
    else c_set_prec_type rawptr solver prec_type true

  let default = function
    | Some x -> x
    | None -> 0

  external c_spbcgs : int -> ('d, 'k) Nvector.t -> ('m, 'nd, 'nk) cptr
    = "sunml_lsolver_spbcgs"

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
    LS {
      rawptr = cptr;
      solver = Spbcgs;
      matrix = None;
      compat = compat;
      check_prec_type = (fun _ -> true);
      attached = false;
    }

  external c_spfgmr : int -> ('d, 'k) Nvector.t -> ('m, 'nd, 'nk) cptr
    = "sunml_lsolver_spfgmr"

  let spfgmr ?maxl ?max_restarts ?gs_type nvec =
    let maxl = default maxl in
    let cptr = c_spfgmr maxl nvec in
    let compat =
      if in_compat_mode
      then begin
          let r = { info with maxl = maxl;
                              gs_type = gs_type;
                              max_restarts = max_restarts;
          } in
          r.set_gs_type <- (fun t -> r.gs_type <- Some t);
          r.set_max_restarts <- (fun t -> r.max_restarts <- Some t);
          r
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
    LS {
      rawptr = cptr;
      solver = Spfgmr;
      matrix = None;
      compat = compat;
      check_prec_type = (fun _ -> true);
      attached = false;
    }

  external c_spgmr : int -> ('d, 'k) Nvector.t -> ('m, 'nd, 'nk) cptr
    = "sunml_lsolver_spgmr"

  let spgmr ?maxl ?max_restarts ?gs_type nvec =
    let maxl = default maxl in
    let cptr = c_spgmr maxl nvec in
    let compat =
      if in_compat_mode
      then begin
          let r = { info with maxl = maxl;
                              gs_type = gs_type;
                              max_restarts = max_restarts;
          } in
          r.set_gs_type <- (fun t -> r.gs_type <- Some t);
          r.set_max_restarts <- (fun t -> r.max_restarts <- Some t);
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
    LS {
      rawptr = cptr;
      solver = Spgmr;
      matrix = None;
      compat = compat;
      check_prec_type = (fun _ -> true);
      attached = false;
    }

  external c_sptfqmr : int -> ('d, 'k) Nvector.t -> ('m, 'nd, 'nk) cptr
    = "sunml_lsolver_sptfqmr"

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
    LS {
      rawptr = cptr;
      solver = Sptfqmr;
      matrix = None;
      compat = compat;
      check_prec_type = (fun _ -> true);
      attached = false;
    }

  external c_pcg : int -> ('d, 'k) Nvector.t -> ('m, 'nd, 'nk) cptr
    = "sunml_lsolver_pcg"

  let pcg ?maxl nvec =
    let maxl = default maxl in
    let cptr = c_pcg maxl nvec in
    let compat =
      if in_compat_mode
      then let r = { info with maxl = maxl } in
           r.set_maxl <- (fun m -> r.maxl <- m); r
      else info
    in
    LS {
      rawptr = cptr;
      solver = Pcg;
      matrix = None;
      compat = compat;
      check_prec_type = (fun _ -> true);
      attached = false;
    }

  module Algorithms = struct (* {{{ *)

    external qr_fact : int
                       -> RealArray2.t
                       -> RealArray.t
                       -> bool
                       -> unit
      = "sunml_spils_qr_fact"

    external qr_sol : int
                      -> RealArray2.t
                      -> RealArray.t
                      -> RealArray.t
                      -> unit
      = "sunml_spils_qr_sol"

    external modified_gs : (('a, 'k) Nvector.t) array
                           -> RealArray2.t
                           -> int
                           -> int
                           -> float
      = "sunml_spils_modified_gs"

    external classical_gs' : (('a, 'k) Nvector.t) array
                             * RealArray2.t
                             * int
                             * int
                             * RealArray.t
                             * (('a, 'k) Nvector.t) array
                             -> float
      = "sunml_spils_classical_gs"

    let classical_gs v h k p s temp =
      if Sundials_configuration.safe then begin
        if k < 1 then
          invalid_arg "classical_gs: k is too small.";
        if Array.length v < k then
          invalid_arg "classical_gs: v is too small (< k).";
        if Array.length temp < k + 1 then
          invalid_arg "classical_gs: temp is too small (< k + 1).";
        Array.iter (Nvector.check v.(0)) v;
        Array.iter (Nvector.check v.(0)) temp
      end;
      classical_gs' (v, h, k, p, s, temp)

  end (* }}} *)
end (* }}} *)

module Custom = struct (* {{{ *)

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

  type ('matrix, 'data, 'kind, 'lsolver) ops = {

      init : 'lsolver -> unit;

      setup : 'lsolver -> 'matrix -> unit;

      solve : 'lsolver -> 'matrix -> 'data -> 'data -> float -> unit;

      set_atimes
      : ('lsolver -> ('data, 'kind) atimesfn -> unit) option;

      set_preconditioner
      : ('lsolver 
         -> psetupfn option
         -> ('data, 'kind) psolvefn option
         -> unit) option;

      set_scaling_vectors
      : ('lsolver -> 'data option -> 'data option -> unit) option;

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
       fun fd -> fset' (LSI.Custom.call_atimes fd)

  let wrap_set_preconditioner fseto ldata =
    match fseto with
    | None -> fun _  ->
              failwith ("internal error: Iterative.Custom.set_preconditioner")
    | Some fset ->
       let fset' = fset ldata in
       fun fd has_setup has_solve ->
       fset'
         (if has_setup
          then Some (fun () -> LSI.Custom.call_psetup fd) else None)
         (if has_solve
          then Some (LSI.Custom.call_psolve fd) else None)

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
             get_work_space = fget_work_space } ldata mat =
    (match Config.sundials_version with
     | 2,_,_ -> raise Config.NotImplementedBySundialsVersion;
     | _ -> ());
    let ops = LSI.Custom.({
        init = (fun () -> finit ldata);
        setup = fsetup ldata;
        solve = fsolve ldata;
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
    let only_ops = LSI.Custom.({
          has_set_atimes          = fset_atimes <> None;
          has_set_preconditioner  = fset_preconditioner <> None;
          has_set_scaling_vectors = fset_scaling_vectors <> None;
          has_get_num_iters       = fget_num_iters <> None;
          has_get_res_norm        = fget_res_norm <> None;
          has_get_res_id          = fget_res_id <> None;
          has_get_work_space      = fget_work_space <> None;
       })
    in
    LS {
       rawptr = c_make_custom 1 ops only_ops;
       solver = Custom (ldata, ops);
       matrix = mat;
       compat = LSI.Iterative.info;
       check_prec_type = (fun _ -> true);
       attached = false;
     }

  let unwrap (LS { solver } :
        ('m, 'data, 'kind, [>`Custom of 'lsolver]) linear_solver) =
    match solver with
    | Custom (ldata, _) -> ldata
    | _ -> assert false (* Guaranteed by type constraint but not
                           inferred by the type system. *)

  type ('matrix, 'data, 'kind, 'lsolver) dls_ops = {
      init : 'lsolver -> unit;

      setup : 'lsolver -> 'matrix -> unit;

      solve : 'lsolver -> 'matrix -> 'data -> 'data -> float -> unit;

      space : ('lsolver -> int * int) option;
    }

  let make_dls { init = fi; setup = fs0; solve = fs; space = fgws}
               ldata mat =
    (match Config.sundials_version with
     | 2,_,_ -> raise Config.NotImplementedBySundialsVersion;
     | _ -> ());
    let ops = LSI.Custom.({
        init = (fun () -> fi ldata);
        setup = fs0 ldata;
        solve = fs ldata;
        set_atimes = (fun _ ->
          failwith "internal error: Direct.Custom.set_atimes");
        set_preconditioner = (fun _ _ ->
          failwith "internal error: Direct.Custom.set_preconditioner");
        set_scaling_vectors = (fun _ _ ->
          failwith "internal error: Direct.Custom.set_scaling_vectors");
        get_num_iters = (fun _ ->
          failwith "internal error: Direct.Custom.get_num_iters");
        get_res_norm = (fun _ ->
          failwith "internal error: Direct.Custom.get_res_norm");
        get_res_id = (fun _ ->
          failwith "internal error: Direct.Custom.get_res_id");
        get_work_space = match fgws with
          | Some f -> (fun _ -> f ldata)
          | None -> (fun () -> (0, 0))
      }) in
    let only_ops = LSI.Custom.({
          has_set_atimes          = false;
          has_set_preconditioner  = false;
          has_set_scaling_vectors = false;
          has_get_num_iters       = false;
          has_get_res_norm        = false;
          has_get_res_id          = false;
          has_get_work_space      = fgws <> None;
       })
    in
    LS {
     rawptr = c_make_custom 0 ops only_ops;
     solver = Custom (ldata, ops);
     matrix = Some mat;
     compat = LSI.Iterative.info;
     check_prec_type = (fun _ -> true);
     attached = false;
   }

end (* }}} *)

(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "sunml_lsolver_init_module"

let _ =
  c_init_module
    (* Exceptions must be listed in the same order as
       lsolver_exn_index.  *)
    [|InvalidLinearSolver;
      UnrecoverableFailure false;
      MatrixNotSquare;
      MatrixVectorMismatch;
      InsufficientStorageUpperBandwidth;
      Invalid_argument ""; (* Standard OCaml exception *)
      ATimesFailure false;
      PSetFailure false;
      PSolveFailure false;
      GSFailure;
      QRSolFailure;
      VectorOpError;
      ResReduced;
      ConvFailure;
      QRfactFailure;
      LUfactFailure;
      PackageFailure false;
      IllegalPrecType;
      InternalFailure ("", 0);
      ZeroInDiagonal 0;
    |]

