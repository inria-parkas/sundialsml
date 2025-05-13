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

open LSI

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

type ('mat, 'data, 'kind, 'tag) t = ('mat, 'data, 'kind, 'tag) linear_solver

(** Alias for linear solvers that are restricted to serial nvectors. *)
type ('mat, 'kind, 'tag) serial_t
  = ('mat, Nvector_serial.data, [>Nvector_serial.kind] as 'kind, 'tag)
      linear_solver

(* Must correspond with linearSolver_ml.h:lsolver_linear_solver_type *)
type linear_solver_type = LSI.linear_solver_type =
  | Direct
  | Iterative
  | MatrixIterative
  | MatrixEmbedded

type linear_solver_id = LSI.linear_solver_id =
  | Band
  | Dense
  | Klu
  | LapackBand
  | LapackDense
  | Pcg
  | Spbcgs
  | Spfgmr
  | Spgmr
  | Sptfqmr
  | Superludist
  | Superlumt
  | CuSolverSp_batchQR
  | MagmaDense
  | OneMKLDense
  | Gingko
  | KokkosDense
  | Custom

exception LinearSolverInUse = LSI.LinearSolverInUse

type 'd atimesfn = 'd LSI.atimesfn
type psetupfn = LSI.psetupfn
type 'd psolvefn = 'd LSI.psolvefn

module Direct = struct (* {{{ *)

  external c_dense
           : 'k Nvector.serial
             -> 'k Matrix.dense
             -> Sundials.Context.t
             -> (Matrix.Dense.t, Nvector_serial.data, 'k) cptr
    = "sunml_lsolver_dense"

  let dense ?context nvec mat =
    let ctx = Sundials_impl.Context.get context in
    LS {
      rawptr = c_dense nvec mat ctx;
      solver = Dense;
      matrix = Some mat;
      compat = LSI.Iterative.info;
      context = ctx;
      check_prec_type = (fun _ -> true);
      ocaml_callbacks = empty_ocaml_callbacks ();
      info_file = None;
      attached = false;
    }

  external c_lapack_dense
           : 'k Nvector.serial
             -> 'k Matrix.dense
             -> Sundials.Context.t
             -> (Matrix.Dense.t, Nvector_serial.data, 'k) cptr
    = "sunml_lsolver_lapack_dense"

  let lapack_dense ?context nvec mat =
    if not Config.lapack_enabled
    then raise Config.NotImplementedBySundialsVersion;
    let ctx = Sundials_impl.Context.get context in
    LS {
      rawptr = c_lapack_dense nvec mat ctx;
      solver = LapackDense;
      matrix = Some mat;
      compat = LSI.Iterative.info;
      context = ctx;
      check_prec_type = (fun _ -> true);
      ocaml_callbacks = empty_ocaml_callbacks ();
      info_file = None;
      attached = false;
    }

  external c_band
           : 'k Nvector.serial
             -> 'k Matrix.band
             -> Sundials.Context.t
             -> (Matrix.Band.t, Nvector_serial.data, 'k) cptr
    = "sunml_lsolver_band"

  let band ?context nvec mat =
    let ctx = Sundials_impl.Context.get context in
    LS {
      rawptr = c_band nvec mat ctx;
      solver = Band;
      matrix = Some mat;
      compat = LSI.Iterative.info;
      context = ctx;
      check_prec_type = (fun _ -> true);
      ocaml_callbacks = empty_ocaml_callbacks ();
      info_file = None;
      attached = false;
    }

  external c_lapack_band
           : 'k Nvector.serial
             -> 'k Matrix.band
             -> Sundials.Context.t
             -> (Matrix.Band.t, Nvector_serial.data, 'k) cptr
    = "sunml_lsolver_lapack_band"

  let lapack_band ?context nvec mat =
    if not Config.lapack_enabled
    then raise Config.NotImplementedBySundialsVersion;
    let ctx = Sundials_impl.Context.get context in
    LS {
      rawptr = c_lapack_band nvec mat ctx;
      solver = LapackBand;
      matrix = Some mat;
      compat = LSI.Iterative.info;
      context = ctx;
      check_prec_type = (fun _ -> true);
      ocaml_callbacks = empty_ocaml_callbacks ();
      info_file = None;
      attached = false;
    }

  module Klu = struct (* {{{ *)
    include LSI.Klu

    external c_klu
             : 'k Nvector.serial
               -> ('s, 'k) Matrix.sparse
               -> Sundials.Context.t
               -> ('s Matrix.Sparse.t, Nvector_serial.data, 'k) cptr
      = "sunml_lsolver_klu"

    let make ?context ?ordering nvec mat =
      if not Config.klu_enabled
      then raise Config.NotImplementedBySundialsVersion;
      let ctx = Sundials_impl.Context.get context in
      let cptr = c_klu nvec mat ctx in
      let info =
        if Sundials_impl.Version.in_compat_mode2
        then let r = { (info()) with ordering = ordering } in
             r.set_ordering <- (fun o -> r.ordering <- Some o); r
        else info ()
      in
      LS {
        rawptr = cptr;
        solver = Klu info;
        matrix = Some mat;
        compat = LSI.Iterative.info;
        context = ctx;
        check_prec_type = (fun _ -> true);
        ocaml_callbacks = empty_ocaml_callbacks ();
        info_file = None;
        attached = false;
      }

    external c_reinit
             : ('s Matrix.Sparse.t, Nvector_serial.data, 'k) cptr
               -> ('s, 'k) Matrix.sparse
               -> unit
      = "sunml_lsolver_klu_reinit"

    let reinit (LS { rawptr = cptr; solver }) mat ?nnz () =
      if Sundials_impl.Version.in_compat_mode2 then
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
      if Sundials_impl.Version.in_compat_mode2 then
        match solver with
        | Klu { set_ordering = f } -> f ordering
        | _ -> assert false
      else c_set_ordering cptr ordering

  end (* }}} *)

  let klu = Klu.make

  module Superlumt = struct (* {{{ *)
    include LSI.Superlumt

    external c_superlumt
             : 'k Nvector.serial
               -> ('s, 'k) Matrix.sparse
               -> int
               -> Sundials.Context.t
               -> ('s Matrix.Sparse.t, Nvector_serial.data, 'k) cptr
      = "sunml_lsolver_superlumt"

    let make ?context ?ordering ~nthreads nvec mat =
      if not Config.superlumt_enabled
         || (Sundials_impl.Version.in_compat_mode2
             && not Matrix.(Sparse.is_csc (unwrap mat)))
      then raise Config.NotImplementedBySundialsVersion;
      let ctx = Sundials_impl.Context.get context in
      let cptr = c_superlumt nvec mat nthreads ctx in
      let info =
        if Sundials_impl.Version.in_compat_mode2
        then let r = { (info nthreads) with ordering = ordering } in
             r.set_ordering <- (fun o -> r.ordering <- Some o); r
        else info nthreads
      in
      LS {
        rawptr = cptr;
        solver = Superlumt info;
        matrix = Some mat;
        compat = LSI.Iterative.info;
        context = ctx;
        check_prec_type = (fun _ -> true);
        ocaml_callbacks = empty_ocaml_callbacks ();
        info_file = None;
        attached = false;
      }

    external c_set_ordering
             : ('s Matrix.Sparse.t, Nvector_serial.data, 'k) cptr
               -> ordering
               -> unit
      = "sunml_lsolver_superlumt_set_ordering"

    let set_ordering (LS { rawptr = cptr; solver }) ordering =
      if Sundials_impl.Version.in_compat_mode2 then
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
    if Sundials_impl.Version.in_compat_mode2 then compat.set_maxl maxl
    else c_set_maxl rawptr solver maxl

  let set_gs_type (LS { rawptr; solver; compat }) gs_type =
    if Sundials_impl.Version.in_compat_mode2 then compat.set_gs_type gs_type
    else c_set_gs_type rawptr solver gs_type

  let set_max_restarts (LS { rawptr; solver; compat }) max_restarts =
    if Sundials_impl.Version.in_compat_mode2 then compat.set_max_restarts max_restarts
    else c_set_max_restarts rawptr solver max_restarts

  let set_prec_type lsolver prec_type =
    let LS { rawptr; solver; compat; check_prec_type } = lsolver in
    if not (check_prec_type prec_type) then raise IllegalPrecType;
    if Sundials_impl.Version.in_compat_mode2 then compat.set_prec_type prec_type
    else impl_set_prec_type rawptr solver prec_type true

  external c_set_print_level
   : ('m, 'nd, 'nk) cptr
     -> ('m, 'nd, 'nk, [> `Iter]) solver_data
     -> int
     -> unit
   = "sunml_lsolver_set_print_level"

  let set_print_level (LS { rawptr; solver; _ }) level =
    c_set_print_level rawptr solver (if level then 1 else 0)

  external c_set_info_file
   : ('m, 'nd, 'nk) cptr
     -> ('m, 'nd, 'nk, [> `Iter]) solver_data
     -> Logfile.t
     -> unit
   = "sunml_lsolver_set_info_file"

  let set_info_file (LS ({ rawptr; solver; _ } as lsdata)) ?print_level file =
    lsdata.info_file <- Some file;
    c_set_info_file rawptr solver file;
    (match print_level with None -> ()
     | Some level -> c_set_print_level rawptr solver (if level then 1 else 0))

  let default = function
    | Some x -> x
    | None -> 0

  external c_spbcgs
    : int -> ('d, 'k) Nvector.t -> Sundials.Context.t -> ('m, 'nd, 'nk) cptr
    = "sunml_lsolver_spbcgs"

  let spbcgs ?context ?maxl nvec =
    let maxl = default maxl in
    let ctx = Sundials_impl.Context.get context in
    let cptr = c_spbcgs maxl nvec ctx in
    let compat =
      if Sundials_impl.Version.in_compat_mode2
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
      context = ctx;
      check_prec_type = (fun _ -> true);
      ocaml_callbacks = empty_ocaml_callbacks ();
      info_file = None;
      attached = false;
    }

  external c_spfgmr
    : int -> ('d, 'k) Nvector.t -> Sundials.Context.t -> ('m, 'nd, 'nk) cptr
    = "sunml_lsolver_spfgmr"

  let spfgmr ?context ?maxl ?max_restarts ?gs_type nvec =
    let maxl = default maxl in
    let ctx = Sundials_impl.Context.get context in
    let cptr = c_spfgmr maxl nvec ctx in
    let compat =
      if Sundials_impl.Version.in_compat_mode2
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
      context = ctx;
      check_prec_type = (fun _ -> true);
      ocaml_callbacks = empty_ocaml_callbacks ();
      info_file = None;
      attached = false;
    }

  external c_spgmr
    : int -> ('d, 'k) Nvector.t -> Sundials.Context.t -> ('m, 'nd, 'nk) cptr
    = "sunml_lsolver_spgmr"

  let spgmr ?context ?maxl ?max_restarts ?gs_type nvec =
    let maxl = default maxl in
    let ctx = Sundials_impl.Context.get context in
    let cptr = c_spgmr maxl nvec ctx in
    let compat =
      if Sundials_impl.Version.in_compat_mode2
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
      context = ctx;
      check_prec_type = (fun _ -> true);
      ocaml_callbacks = empty_ocaml_callbacks ();
      info_file = None;
      attached = false;
    }

  external c_sptfqmr
    : int -> ('d, 'k) Nvector.t -> Sundials.Context.t -> ('m, 'nd, 'nk) cptr
    = "sunml_lsolver_sptfqmr"

  let sptfqmr ?context ?maxl nvec =
    let maxl = default maxl in
    let ctx = Sundials_impl.Context.get context in
    let cptr = c_sptfqmr maxl nvec ctx in
    let compat =
      if Sundials_impl.Version.in_compat_mode2
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
      context = ctx;
      check_prec_type = (fun _ -> true);
      ocaml_callbacks = empty_ocaml_callbacks ();
      info_file = None;
      attached = false;
    }

  external c_pcg
    : int -> ('d, 'k) Nvector.t -> Sundials.Context.t -> ('m, 'nd, 'nk) cptr
    = "sunml_lsolver_pcg"

  let pcg ?context ?maxl nvec =
    let maxl = default maxl in
    let ctx = Sundials_impl.Context.get context in
    let cptr = c_pcg maxl nvec ctx in
    let compat =
      if Sundials_impl.Version.in_compat_mode2
      then let r = { info with maxl = maxl } in
           r.set_maxl <- (fun m -> r.maxl <- m); r
      else info
    in
    LS {
      rawptr = cptr;
      solver = Pcg;
      matrix = None;
      compat = compat;
      context = ctx;
      check_prec_type = (fun _ -> true);
      ocaml_callbacks = empty_ocaml_callbacks ();
      info_file = None;
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

  type ('d, 'k) atimesfn =
    ('d, 'k) Nvector.t -> ('d, 'k) Nvector.t -> unit

  type psetupfn = unit -> unit

  type ('d, 'k) psolvefn =
    ('d, 'k) Nvector.t -> ('d, 'k) Nvector.t -> float -> bool -> unit

  type ('matrix, 'data, 'kind, 'lsolver) ops = {
      solver_type : linear_solver_type;
      solver_id   : linear_solver_id;

      init : ('lsolver -> unit) option;

      setup : ('lsolver -> 'matrix option -> unit) option;

      solve : 'lsolver -> 'matrix option -> 'data -> 'data -> float -> unit;

      set_atimes
      : ('lsolver -> ('data, 'kind) atimesfn -> unit) option;

      set_preconditioner
      : ('lsolver 
         -> psetupfn option
         -> ('data, 'kind) psolvefn option
         -> unit) option;

      set_scaling_vectors
      : ('lsolver -> 'data option -> 'data option -> unit) option;

      set_zero_guess : ('lsolver -> bool -> unit) option;

      get_num_iters : ('lsolver -> int) option;

      get_res_norm : ('lsolver -> float) option;

      get_res_id : ('lsolver -> ('data, 'kind) Nvector.t) option;

      get_last_flag : ('lsolver -> int) option;

      get_work_space : ('lsolver -> int * int) option;

      set_prec_type
      : ('lsolver -> Iterative.preconditioning_type -> unit) option;
    }

  let make_ops ?(solver_id=Custom) ?init ?setup
               ?set_atimes ?set_preconditioner ?set_scaling_vectors
               ?set_zero_guess ?get_num_iters ?get_res_norm ?get_res_id
               ?get_last_flag ?get_work_space ?set_prec_type
               ~solver_type ~solve () =
    { solver_type; solver_id  ; init; setup; solve; set_atimes;
      set_preconditioner; set_scaling_vectors; set_zero_guess;
      get_num_iters; get_res_norm; get_res_id; get_last_flag;
      get_work_space; set_prec_type; }

  let wrap_set_atimes fseto =
    match fseto with
    | None ->
       fun _ _  -> failwith ("internal error: Iterative.Custom.set_atimes")
    | Some fset ->
       fun ldata fd -> fset ldata (LSI.Custom.call_atimes fd)

  let wrap_set_preconditioner fseto =
    match fseto with
    | None -> fun _ _  ->
              failwith ("internal error: Iterative.Custom.set_preconditioner")
    | Some fset ->
       fun ldata fd has_setup has_solve ->
       fset ldata
              (if has_setup
               then Some (fun () -> LSI.Custom.call_psetup fd) else None)
              (if has_solve
               then Some (LSI.Custom.call_psolve fd) else None)

  let mapignore fo =
    match fo with
    | None -> (fun _ _ -> ())
    | Some f -> f

  let mapo s fo =
    match fo with
    | None -> (fun _ -> failwith ("internal error: Iterative.Custom." ^ s))
    | Some f -> f

  let weak_wrap x =
    let wx = Weak.create 1 in
    Weak.set wx 0 (Some x);
    wx

  [@@@warning "-45"]
  let make { solver_type = stype;
             solver_id = sid;
             init = finit;
             setup = fsetup;
             solve = fsolve;
             set_atimes = fset_atimes;
             set_preconditioner = fset_preconditioner;
             set_scaling_vectors = fset_scaling_vectors;
             set_zero_guess = fset_zero_guess;
             get_num_iters = fget_num_iters;
             get_res_norm = fget_res_norm;
             get_res_id = fget_res_id;
             get_last_flag = fget_last_flag;
             get_work_space = fget_work_space;
             set_prec_type = fset_prec_type } ?context ldata omat =
    (match Config.sundials_version with
     | 2,_,_ -> raise Config.NotImplementedBySundialsVersion;
     | _ -> ());
    let ops = LSI.Custom.({
        init = (match finit with
                | None -> (fun _ -> ())
                | Some f -> f);
        setup = (match fsetup with
                 | None -> (fun _ _ -> ())
                 | Some f -> f);
        solve = fsolve;
        set_atimes = wrap_set_atimes fset_atimes;
        set_preconditioner =
          wrap_set_preconditioner fset_preconditioner;
        set_scaling_vectors =
          mapo "set_scaling_vectors" fset_scaling_vectors;
        set_zero_guess = mapo "set_zero_guess" fset_zero_guess;
        get_id = (fun _ -> sid);
        get_num_iters = mapo "get_num_iters" fget_num_iters;
        get_res_norm = mapo "get_res_norm" fget_res_norm;
        get_res_id = mapo "get_res_id" fget_res_id;
        get_last_flag = mapo "get_last_flag" fget_last_flag;
        get_work_space = mapo "get_work_space" fget_work_space;
        (* set_prec_type is only every called from OCaml and never from C. *)
        set_prec_type = mapignore fset_prec_type;
      }) in
    let only_ops = LSI.Custom.({
          has_init                = finit <> None;
          has_setup               = fsetup <> None
                                    && stype <> MatrixEmbedded;
          has_set_atimes          = fset_atimes <> None
                                    && stype <> MatrixEmbedded;
          has_set_preconditioner  = fset_preconditioner <> None;
          has_set_scaling_vectors = fset_scaling_vectors <> None;
          has_set_zero_guess      = fset_zero_guess <> None;
          has_get_num_iters       = fget_num_iters <> None;
          has_get_res_norm        = fget_res_norm <> None;
          has_get_res_id          = fget_res_id <> None;
          has_get_last_flag       = fget_last_flag <> None;
          has_get_work_space      = fget_work_space <> None;
       })
    in
    let ldata_and_ops = (ldata, ops) in
    let ctx = Sundials_impl.Context.get context in
    (LS {
       rawptr = c_make_custom stype (weak_wrap ldata_and_ops) only_ops ctx;
       solver = Custom ldata_and_ops;
       matrix = omat;
       compat = LSI.Iterative.info;
       context = ctx;
       check_prec_type = (fun _ -> true);
       ocaml_callbacks = empty_ocaml_callbacks ();
       info_file = None;
       attached = false;
     })
  [@@@warning "+45"]

  let make_with_matrix ({ solver_type; _ } as ops) ?context ldata (mat : 'matrix) =
    if solver_type = MatrixEmbedded
      then invalid_arg "invalid solver_type when matrix is given";
    make ops ?context ldata (Some mat)

  let make_without_matrix ({ solver_type; _ } as ops) ?context ldata =
    if solver_type = Direct
      then invalid_arg "invalid solver_type when matrix is not given";
    if Sundials_impl.Version.lt580 && solver_type = MatrixEmbedded
      then raise Config.NotImplementedBySundialsVersion;
    make ops ?context ldata None

  let unwrap (LS { solver } :
        ('m, 'data, 'kind, [>`Custom of 'lsolver]) linear_solver) =
    match solver with
    | Custom (ldata, _) -> ldata
    | _ -> assert false (* Guaranteed by type constraint but not
                           inferred by the type system. *)

  type ('matrix, 'data, 'kind, 'lsolver) dls_ops = {
      init : ('lsolver -> unit) option;

      setup : ('lsolver -> 'matrix -> unit) option;

      solve : 'lsolver -> 'matrix -> 'data -> 'data -> float -> unit;

      space : ('lsolver -> int * int) option;
    }

  [@@@warning "-45"]
  let make_dls { init = fi; setup = fs0; solve = fsolve; space = fgws} ?context
               ldata mat =
    (match Config.sundials_version with
     | 2,_,_ -> raise Config.NotImplementedBySundialsVersion;
     | _ -> ());
    let ops = LSI.Custom.({
        init = (match fi with
                | None -> (fun _ -> ())
                | Some f -> f);
        setup = (match fs0 with
                 | None -> (fun _ _ -> ())
                 | Some f -> (fun ls omat -> f ls (Option.get omat)));
        solve = (fun ls omat v1 v2 -> fsolve ls (Option.get omat) v1 v2);
        set_atimes = (fun _ _ ->
          failwith "internal error: Direct.Custom.set_atimes");
        set_preconditioner = (fun _ _ ->
          failwith "internal error: Direct.Custom.set_preconditioner");
        set_scaling_vectors = (fun _ _ ->
          failwith "internal error: Direct.Custom.set_scaling_vectors");
        set_zero_guess = (fun _ ->
          failwith "internal error: Direct.Custom.set_zero_guess");
        get_id = (fun _ -> Custom);
        get_num_iters = (fun _ ->
          failwith "internal error: Direct.Custom.get_num_iters");
        get_res_norm = (fun _ ->
          failwith "internal error: Direct.Custom.get_res_norm");
        get_res_id = (fun _ ->
          failwith "internal error: Direct.Custom.get_res_id");
        get_last_flag = (fun _ ->
          failwith "internal error: Direct.Custom.get_last_flag");
        get_work_space = (match fgws with
          | Some f -> f
          | None -> (fun _ -> (0, 0)));
        set_prec_type = (fun _ _ -> ());
      }) in
    let only_ops = LSI.Custom.({
          has_init                = fi <> None;
          has_setup               = fs0 <> None;
          has_set_atimes          = false;
          has_set_preconditioner  = false;
          has_set_scaling_vectors = false;
          has_set_zero_guess      = false;
          has_get_num_iters       = false;
          has_get_res_norm        = false;
          has_get_res_id          = false;
          has_get_last_flag       = false;
          has_get_work_space      = fgws <> None;
       })
    in
    let ldata_and_ops = (ldata, ops) in
    let ctx = Sundials_impl.Context.get context in
    LS {
     rawptr = c_make_custom Direct (weak_wrap ldata_and_ops) only_ops ctx;
     solver = Custom ldata_and_ops;
     matrix = Some mat;
     compat = LSI.Iterative.info;
     context = ctx;
     check_prec_type = (fun _ -> true);
     ocaml_callbacks = empty_ocaml_callbacks ();
     info_file = None;
     attached = false;
   }
  [@@@warning "+45"]

end (* }}} *)

external c_get_type : ('m, 'd, 'k) cptr -> linear_solver_type
  = "sunml_lsolver_get_type"

let get_type (LS { rawptr }) = c_get_type rawptr

external c_get_id : ('m, 'd, 'k) cptr -> linear_solver_id
  = "sunml_lsolver_get_id"

let get_id (LS { rawptr }) = c_get_id rawptr

external c_set_atimes
  : ('m, 'd, 'k) cptr
    -> ('d, 'k) ocaml_callbacks Sundials_impl.Vptr.vptr
    -> unit
  = "sunml_lsolver_set_atimes"

let set_atimes (LS { rawptr; ocaml_callbacks }) fn =
  (Sundials_impl.Vptr.unwrap ocaml_callbacks).ocaml_atimes <- fn;
  c_set_atimes rawptr ocaml_callbacks

external c_set_preconditioner
  : ('m, 'd, 'k) cptr
    -> ('d, 'k) ocaml_callbacks Sundials_impl.Vptr.vptr
    -> unit
  = "sunml_lsolver_set_preconditioner"

let set_preconditioner (LS { rawptr; ocaml_callbacks }) psetupf psolvef =
  let cb = Sundials_impl.Vptr.unwrap ocaml_callbacks in
  cb.ocaml_psetup <- psetupf;
  cb.ocaml_psolve <- psolvef;
  c_set_preconditioner rawptr ocaml_callbacks

external c_set_scaling_vectors
  : ('m, 'd, 'k) cptr
    -> ('d, 'k) Nvector.t
    -> ('d, 'k) Nvector.t
    -> unit
  = "sunml_lsolver_set_scaling_vectors"

let set_scaling_vectors (LS { rawptr; ocaml_callbacks }) s1 s2 =
  let cb = Sundials_impl.Vptr.unwrap ocaml_callbacks in
  cb.scaling_vector1 <- Some s1; (* prevent garbage collection *)
  cb.scaling_vector2 <- Some s2; (* prevent garbage collection *)
  c_set_scaling_vectors rawptr s1 s2

external c_set_zero_guess
  : ('m, 'd, 'k) cptr -> bool -> unit
  = "sunml_lsolver_set_zero_guess"

let set_zero_guess (LS { rawptr; _ }) onoff =
  c_set_zero_guess rawptr onoff

external c_initialize : ('m, 'd, 'k) cptr -> unit
  = "sunml_lsolver_initialize"

let init (LS { rawptr }) = c_initialize rawptr

external c_setup :
     ('m, 'nd, 'nk) cptr
  -> ('a, 'm, 'd, 'k) Matrix.t option
  -> unit
  = "sunml_lsolver_setup"

let setup (LS { rawptr }) = c_setup rawptr

external c_solve :
     ('m, 'd, 'k) cptr
  -> ('a, 'm, 'd, 'k) Matrix.t option
  -> ('d, 'k) Nvector.t
  -> ('d, 'k) Nvector.t
  -> float
  -> unit
  = "sunml_lsolver_solve"

let solve (LS { rawptr }) = c_solve rawptr

external c_iters : ('m, 'd, 'k) cptr -> int
  = "sunml_lsolver_iters"

let get_num_iters (LS { rawptr }) = c_iters rawptr

external c_res_norm : ('m, 'd, 'k) cptr -> float
  = "sunml_lsolver_res_norm"

let get_res_norm (LS { rawptr }) = c_res_norm rawptr

external c_res_id : ('m, 'd, 'k) cptr -> 'd
  = "sunml_lsolver_res_id"

let get_res_id (LS { rawptr }) = c_res_id rawptr

external c_last_flag : ('m, 'd, 'k) cptr -> int
  = "sunml_lsolver_last_flag"

let get_last_flag (LS { rawptr }) = c_last_flag rawptr

external c_space : ('m, 'd, 'k) cptr -> int * int
  = "sunml_lsolver_space"

let get_work_space (LS { rawptr }) = c_space rawptr

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

