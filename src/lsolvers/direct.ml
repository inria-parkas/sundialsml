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

include Lsolver_impl.Direct

type ('m, 'nk, 'tag) serial_t
  = ('m, Nvector_serial.data, [>Nvector_serial.kind] as 'nk, 'tag) t

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
  include Lsolver_impl.Klu

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
      then let r = { info with ordering = ordering } in
           r.set_ordering <- (fun o -> r.ordering <- Some o); r
      else info
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
  include Lsolver_impl.Superlumt

  external c_superlumt
    : 'k Nvector_serial.any
      -> ('s, 'k) Matrix.sparse
      -> int
      -> ('s Matrix.Sparse.t, Nvector_serial.data, 'k) cptr
    = "ml_lsolver_superlumt"

  let make ?ordering ~nthreads nvec mat =
    if not Sundials_config.superlumt_enabled
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

  type 'lsolver tag = 'lsolver Lsolver_impl.Custom.tag

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
    match Sundials.sundials_version with
    | 2,_,_ -> raise Sundials.NotImplementedBySundialsVersion;
    | _ -> ();
    let ops = Lsolver_impl.Custom.({
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
            get_work_space = (fun _ ->
                failwith "internal error: Direct.Custom.get_work_space");
          })
    in
    S {
        rawptr = Lsolver_impl.Direct.(c_make_custom 0 ops only_ops);
        solver = Custom (ldata, ops);
        matrix = mat;
        attached = false;
      }

  let unwrap (S { Lsolver_impl.Direct.solver = Custom (ldata, _) }) = ldata

end (* }}} *)

