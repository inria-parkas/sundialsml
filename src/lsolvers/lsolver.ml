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

exception UnrecoverableFailure of bool
exception MatrixNotSquare

module Direct = struct
  include Lsolver_impl.Direct

  type ('m, 'nk) serial_t
    = ('m, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

  external c_dense : 'k Nvector_serial.any -> 'k Matrix.dense -> cptr
    = "ml_lsolver_dense"

  let dense nvec mat = {
      rawptr = c_dense nvec mat;
      compat = Dense;
    }

  external c_lapack_dense : 'k Nvector_serial.any -> 'k Matrix.dense -> cptr
    = "ml_lsolver_lapack_dense"

  let lapack_dense nvec mat =
    if not Sundials_config.lapack_enabled
      then raise Sundials.NotImplementedBySundialsVersion;
    {
      rawptr = c_lapack_dense nvec mat;
      compat = LapackDense;
    }

  external c_band : 'k Nvector_serial.any -> 'k Matrix.band -> cptr
    = "ml_lsolver_band"

  let band nvec mat = {
      rawptr = c_band nvec mat;
      compat = Band;
    }

  external c_lapack_band : 'k Nvector_serial.any -> 'k Matrix.band -> cptr
    = "ml_lsolver_lapack_band"

  let lapack_band nvec mat =
    if not Sundials_config.lapack_enabled
      then raise Sundials.NotImplementedBySundialsVersion;
    {
      rawptr = c_lapack_band nvec mat;
      compat = LapackBand;
    }

  module Klu = struct
    include Lsolver_impl.Klu

    external c_klu : 'k Nvector_serial.any -> ('s, 'k) Matrix.sparse -> cptr
      = "ml_lsolver_klu"

    let make ?ordering nvec mat =
      if not Sundials_config.klu_enabled
        then raise Sundials.NotImplementedBySundialsVersion;
      let cptr = c_klu nvec mat in
      let compat =
        if in_compat_mode
        then let r = { info with ordering = ordering } in
             r.set_ordering <- (fun o -> r.ordering <- Some o); r
        else info
      in
      { rawptr = cptr;
        compat = Klu compat; }

    external c_reinit : cptr -> ('s, 'k) Matrix.sparse -> unit
      = "ml_lsolver_klu_reinit"

    let reinit { rawptr = cptr; compat } mat ?nnz () =
      if in_compat_mode then
        match compat with
        | Klu { reinit = f } ->
            let smat = Matrix.unwrap mat in
            f (fst (Matrix.Sparse.size smat)) nnz
        | _ -> assert false
      else begin
        (match nnz with
         | Some n -> Matrix.Sparse.upsize n (Matrix.unwrap mat)
         | None -> ());
        c_reinit cptr mat
      end

    external c_set_ordering : cptr -> ordering -> unit
      = "ml_lsolver_klu_set_ordering"

    let set_ordering { rawptr = cptr; compat } ordering =
      if in_compat_mode then
        match compat with
        | Klu { set_ordering = f } -> f ordering
        | _ -> assert false
      else c_set_ordering cptr ordering
  end

  module Superlumt = struct
    include Lsolver_impl.Superlumt

    external c_superlumt
      : 'k Nvector_serial.any -> ('s, 'k) Matrix.sparse -> int -> cptr
      = "ml_lsolver_superlumt"

    let make ?ordering ~nthreads nvec mat =
      if not Sundials_config.superlumt_enabled
        then raise Sundials.NotImplementedBySundialsVersion;
      let cptr = c_superlumt nvec mat nthreads in
      let compat =
        if in_compat_mode
        then let r = { info with ordering = ordering } in
             r.set_ordering <- (fun o -> r.ordering <- Some o); r
        else info
      in
      { rawptr = cptr;
        compat = Superlumt compat; }

    external c_set_ordering : cptr -> ordering -> unit
      = "ml_lsolver_superlumt_set_ordering"

    let set_ordering { rawptr = cptr; compat } ordering =
      if in_compat_mode then
        match compat with
        | Superlumt { set_ordering = f } -> f ordering
        | _ -> assert false
      else c_set_ordering cptr ordering
  end
end

module Iterative = struct
  include Lsolver_impl.Iterative

  exception ATimesFailure of bool
  exception PSetFailure of bool
  exception PSolveFailure of bool
  exception GSFailure
  exception QRSolFailure
  exception ResReduced
  exception ConvFailure
  exception QRfactFailure
  exception LUfactFailure

  external c_set_maxl : cptr -> solver -> int -> unit
    = "ml_lsolver_set_maxl"

  external c_set_gs_type : cptr -> solver -> gramschmidt_type -> unit
    = "ml_lsolver_set_gs_type"

  external c_set_max_restarts : cptr -> solver -> int -> unit
    = "ml_lsolver_set_max_restarts"

  let set_maxl { rawptr = cptr; solver; compat } maxl =
    if in_compat_mode then compat.set_maxl maxl
    else c_set_maxl cptr solver maxl

  let set_gs_type { rawptr = cptr; solver; compat } gs_type =
    if in_compat_mode then compat.set_gs_type gs_type
    else c_set_gs_type cptr solver gs_type

  let set_max_restarts { rawptr = cptr; solver; compat } max_restarts =
    if in_compat_mode then compat.set_max_restarts max_restarts
    else c_set_max_restarts cptr solver max_restarts

  let default = function
    | Some x -> x
    | None -> 0

  external c_spbcgs : int -> ('d, 'k) Nvector.t -> cptr
    = "ml_lsolver_spbcgs"

  let spbcgs ?maxl nvec =
    let cptr = c_spbcgs (default maxl) nvec in
    let compat =
      if in_compat_mode
      then let r = { info with maxl = maxl } in
           r.set_maxl <- (fun m -> r.maxl <- Some m); r
      else info
    in
    { rawptr = cptr;
      solver = Spbcgs;
      compat = compat }

  external c_spfgmr : int -> ('d, 'k) Nvector.t -> cptr
    = "ml_lsolver_spfgmr"

  let spfgmr ?maxl ?max_restarts ?gs_type nvec =
    let cptr = c_spfgmr (default maxl) nvec in
    let compat =
      if in_compat_mode
      then let r = { info with max_restarts = max_restarts;
                               gs_type = gs_type } in
           r.set_max_restarts <- (fun m -> r.max_restarts <- Some m);
           r.set_gs_type      <- (fun t -> r.gs_type <- Some t); r
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
      compat = compat }

  external c_spgmr : int -> ('d, 'k) Nvector.t -> cptr
    = "ml_lsolver_spgmr"

  let spgmr ?maxl ?max_restarts ?gs_type nvec =
    let cptr = c_spgmr (default maxl) nvec in
    let compat =
      if in_compat_mode
      then let r = { info with max_restarts = max_restarts;
                               gs_type      = gs_type } in
           r.set_max_restarts <- (fun m -> r.max_restarts <- Some m);
           r.set_gs_type      <- (fun t -> r.gs_type <- Some t); r
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
      compat = compat }

  external c_sptfqmr : int -> ('d, 'k) Nvector.t -> cptr
    = "ml_lsolver_sptfqmr"

  let sptfqmr ?maxl nvec =
    let cptr = c_sptfqmr (default maxl) nvec in
    let compat =
      if in_compat_mode
      then let r = { info with maxl = maxl } in
           r.set_maxl <- (fun m -> r.maxl <- Some m); r
      else info
    in
    { rawptr = cptr;
      solver = Sptfqmr;
      compat = compat }

  external c_pcg : int -> ('d, 'k) Nvector.t -> cptr
    = "ml_lsolver_pcg"

  let pcg ?maxl nvec =
    let cptr = c_pcg (default maxl) nvec in
    let compat =
      if in_compat_mode
      then let r = { info with maxl = maxl } in
           r.set_maxl <- (fun m -> r.maxl <- Some m); r
      else info
    in
    { rawptr = cptr;
      solver = Pcg;
      compat = compat }
end

(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "ml_lsolver_init_module"

let _ =
  c_init_module
    (* Exceptions must be listed in the same order as
       lsolver_exn_index.  *)
    [|UnrecoverableFailure false;
      MatrixNotSquare;
      Iterative.ATimesFailure false;
      Iterative.PSetFailure false;
      Iterative.PSolveFailure false;
      Iterative.GSFailure;
      Iterative.QRSolFailure;
      Iterative.ResReduced;
      Iterative.ConvFailure;
      Iterative.QRfactFailure;
      Iterative.LUfactFailure;
    |]

