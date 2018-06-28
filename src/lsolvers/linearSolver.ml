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

exception LinearSolverInUse = LinearSolver_impl.LinearSolverInUse
exception UnrecoverableFailure of bool
exception MatrixNotSquare
exception MatrixVectorMismatch
exception InsufficientStorageUpperBandwidth
exception ATimesFailure of bool
exception PSetFailure of bool
exception PSolveFailure of bool
exception GSFailure
exception QRSolFailure
exception ResReduced
exception ConvFailure
exception QRfactFailure
exception LUfactFailure
exception PackageFailure of bool
exception IllegalPrecType
exception InternalFailure of (string * int)

(* "Simulate" Linear Solvers in Sundials < 3.0.0 *)
let in_compat_mode =
  match Sundials.sundials_version with
  | 2,_,_ -> true
  | _ -> false

module Direct = struct (* {{{ *)
  include LinearSolver_impl.Direct

  type ('m, 'nk, 'tag) serial_linear_solver
    = ('m, Nvector_serial.data, [>Nvector_serial.kind] as 'nk, 'tag)
        linear_solver

  external c_dense
           : 'k Nvector_serial.any
             -> 'k Matrix.dense
             -> (Matrix.Dense.t, Nvector_serial.data, 'k) cptr
    = "ml_lsolver_dense"

  let dense nvec mat = S {
                           rawptr = c_dense nvec mat;
                           solver = Dense;
                           matrix = mat;
                           attached = false;
                         }

  external c_lapack_dense
           : 'k Nvector_serial.any
             -> 'k Matrix.dense
             -> (Matrix.Dense.t, Nvector_serial.data, 'k) cptr
    = "ml_lsolver_lapack_dense"

  let lapack_dense nvec mat =
    if not Sundials_config.lapack_enabled
    then raise Sundials.NotImplementedBySundialsVersion;
    S {
        rawptr = c_lapack_dense nvec mat;
        solver = LapackDense;
        matrix = mat;
        attached = false;
      }

  external c_band
           : 'k Nvector_serial.any
             -> 'k Matrix.band
             -> (Matrix.Band.t, Nvector_serial.data, 'k) cptr
    = "ml_lsolver_band"

  let band nvec mat = S {
                          rawptr = c_band nvec mat;
                          solver = Band;
                          matrix = mat;
                          attached = false;
                        }

  external c_lapack_band
           : 'k Nvector_serial.any
             -> 'k Matrix.band
             -> (Matrix.Band.t, Nvector_serial.data, 'k) cptr
    = "ml_lsolver_lapack_band"

  let lapack_band nvec mat =
    if not Sundials_config.lapack_enabled
    then raise Sundials.NotImplementedBySundialsVersion;
    S {
        rawptr = c_lapack_band nvec mat;
        solver = LapackBand;
        matrix = mat;
        attached = false;
      }

  module Klu = struct (* {{{ *)
    include LinearSolver_impl.Klu

    external c_klu
             : 'k Nvector_serial.any
               -> ('s, 'k) Matrix.sparse
               -> ('s Matrix.Sparse.t, Nvector_serial.data, 'k) cptr
      = "ml_lsolver_klu"

    let make ?ordering nvec mat =
      if not Sundials_config.klu_enabled
      then raise Sundials.NotImplementedBySundialsVersion;
      let cptr = c_klu nvec mat in
      let info =
        if in_compat_mode
        then let r = { (info()) with ordering = ordering } in
             r.set_ordering <- (fun o -> r.ordering <- Some o); r
        else info ()
      in
      S { rawptr = cptr;
          solver = Klu info;
          matrix = mat;
          attached = false; }

    external c_reinit
             : ('s Matrix.Sparse.t, Nvector_serial.data, 'k) cptr
               -> ('s, 'k) Matrix.sparse
               -> unit
      = "ml_lsolver_klu_reinit"

    let reinit (S { rawptr = cptr; solver }) mat ?nnz () =
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
      = "ml_lsolver_klu_set_ordering"

    let set_ordering (S { rawptr = cptr; solver }) ordering =
      if in_compat_mode then
        match solver with
        | Klu { set_ordering = f } -> f ordering
        | _ -> assert false
      else c_set_ordering cptr ordering

  end (* }}} *)

  let klu = Klu.make

  module Superlumt = struct (* {{{ *)
    include LinearSolver_impl.Superlumt

    external c_superlumt
             : 'k Nvector_serial.any
               -> ('s, 'k) Matrix.sparse
               -> int
               -> ('s Matrix.Sparse.t, Nvector_serial.data, 'k) cptr
      = "ml_lsolver_superlumt"

    let make ?ordering ~nthreads nvec mat =
      if not Sundials_config.superlumt_enabled
         || (in_compat_mode && not Matrix.(Sparse.is_csc (unwrap mat)))
      then raise Sundials.NotImplementedBySundialsVersion;
      let cptr = c_superlumt nvec mat nthreads in
      let info =
        if in_compat_mode
        then let r = { (info nthreads) with ordering = ordering } in
             r.set_ordering <- (fun o -> r.ordering <- Some o); r
        else info nthreads
      in
      S { rawptr = cptr;
          solver = Superlumt info;
          matrix = mat;
          attached = false; }

    external c_set_ordering
             : ('s Matrix.Sparse.t, Nvector_serial.data, 'k) cptr
               -> ordering
               -> unit
      = "ml_lsolver_superlumt_set_ordering"

    let set_ordering (S { rawptr = cptr; solver }) ordering =
      if in_compat_mode then
        match solver with
        | Superlumt { set_ordering = f } -> f ordering
        | _ -> assert false
      else c_set_ordering cptr ordering

  end (* }}} *)

  let superlumt = Superlumt.make

  module Custom = struct (* {{{ *)

    type 'lsolver tag = 'lsolver LinearSolver_impl.Custom.tag

    type ('matrix, 'data, 'kind, 'lsolver) ops = {
        init : 'lsolver -> unit;

        setup : 'lsolver -> 'matrix -> unit;

        solve : 'lsolver
                -> 'matrix
                -> ('data, 'kind) Nvector.t
                -> ('data, 'kind) Nvector.t
                -> unit;

        get_work_space : ('lsolver -> int * int) option;
      }

    let make { init = fi; setup = fs0; solve = fs; get_work_space = fgws}
             ldata mat =
      (match Sundials.sundials_version with
       | 2,_,_ -> raise Sundials.NotImplementedBySundialsVersion;
       | _ -> ());
      let ops = LinearSolver_impl.Custom.({
                  init = (fun () -> fi ldata);
                  setup = fs0 ldata;
                  solve = (fun a x b tol -> fs ldata a x b);
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
                })
             in
             S {
                 rawptr = LinearSolver_impl.Direct.(c_make_custom 0 ops only_ops);
                 solver = Custom (ldata, ops);
                 matrix = mat;
                 attached = false;
               }

    let unwrap (S { LinearSolver_impl.Direct.solver = Custom (ldata, _) }) =
      ldata

  end (* }}} *)
end (* }}} *)

module Iterative = struct (* {{{ *)
  include LinearSolver_impl.Iterative

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
    if in_compat_mode then compat.set_max_restarts max_restarts
    else c_set_max_restarts rawptr solver max_restarts

  let set_prec_type { rawptr; solver; compat; check_prec_type } prec_type =
    if not (check_prec_type prec_type) then raise IllegalPrecType;
    if in_compat_mode then compat.set_prec_type prec_type
    else c_set_prec_type rawptr solver prec_type true

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
         fun fd -> fset' (LinearSolver_impl.Custom.call_atimes fd)

    let wrap_set_preconditioner fseto ldata =
      match fseto with
      | None -> fun _  ->
                failwith ("internal error: Iterative.Custom.set_preconditioner")
      | Some fset ->
         let fset' = fset ldata in
         fun fd has_setup has_solve ->
         fset'
           (if has_setup
            then Some (fun () -> LinearSolver_impl.Custom.call_psetup fd) else None)
           (if has_solve
            then Some (LinearSolver_impl.Custom.call_psolve fd) else None)

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
      (match Sundials.sundials_version with
       | 2,_,_ -> raise Sundials.NotImplementedBySundialsVersion;
       | _ -> ());
      let ops = LinearSolver_impl.Custom.({
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
      let only_ops = LinearSolver_impl.Custom.({
                        has_set_atimes          = fset_atimes <> None;
                        has_set_preconditioner  = fset_preconditioner <> None;
                        has_set_scaling_vectors = fset_scaling_vectors <> None;
                        has_get_num_iters       = fget_num_iters <> None;
                        has_get_res_norm        = fget_res_norm <> None;
                        has_get_res_id          = fget_res_id <> None;
                        has_get_work_space      = fget_work_space <> None;
                     })
      in LinearSolver_impl.Iterative.({
           rawptr = c_make_custom 1 ops only_ops;
           solver = Custom (ldata, ops);
           compat = info;
           check_prec_type = (fun _ -> true);
           attached = false;
         })

    let unwrap { LinearSolver_impl.Iterative.solver = Custom (ldata, _) } =
      ldata

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
end (* }}} *)

(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "ml_lsolver_init_module"

let _ =
  c_init_module
    (* Exceptions must be listed in the same order as
       lsolver_exn_index.  *)
    [|UnrecoverableFailure false;
      MatrixNotSquare;
      MatrixVectorMismatch;
      InsufficientStorageUpperBandwidth;
      Invalid_argument ""; (* Standard OCaml exception *)
      ATimesFailure false;
      PSetFailure false;
      PSolveFailure false;
      GSFailure;
      QRSolFailure;
      ResReduced;
      ConvFailure;
      QRfactFailure;
      LUfactFailure;
      PackageFailure false;
      IllegalPrecType;
      InternalFailure ("", 0);
    |]

