(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2018 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

open Sundials

exception Invalidated
exception IncompatibleArguments
exception ZeroDiagonalElement of int

let check_valid =
  match Config.sundials_version with
  | 2,_,_ -> Sundials_configuration.safe
  | _ -> false

(* Must correspond with matrix_ml.h:mat_matrix_content_index *)
type [@warning "-69"] ('data, 'cptr) matrix_content = {
  mutable payload : 'data;
  rawptr  : 'cptr;
  mutable valid : bool;
}

(* Must correspond with matrix_ml.h:mat_matrix_ops_index *)
type ('m, 'd) matrix_ops = {
  m_clone        : 'm -> 'm;

  m_zero         : 'm -> unit;

  m_copy         : 'm -> 'm -> unit;

  m_scale_add    : float -> 'm -> 'm -> unit;

  m_scale_addi   : float -> 'm -> unit;

  m_matvec_setup : ('m -> unit) option;

  m_matvec       : 'm -> 'd -> 'd -> unit;

  m_space        : 'm -> int * int;
}

module Dense = struct (* {{{ *)

  type data = RealArray2.data
  type cptr

  type t = (data, cptr) matrix_content

  external c_create : int -> int -> t
    = "sunml_matrix_dense_create"

  let create i j =
    if Sundials_configuration.safe then begin
      if i <= 0 then invalid_arg "number of rows";
      if j <= 0 then invalid_arg "number of cols"
    end;
    c_create i j

  let unwrap { payload } = payload

  let make m n v =
    let { payload } as r = create m n in
    Bigarray.Array2.fill payload v;
    r

  let invalidate v =
    if check_valid then v.valid <- false

  let size { payload; valid } =
    if check_valid && not valid then raise Invalidated;
    Bigarray.Array2.(dim2 payload, dim1 payload)

  let pp fmt { payload = d; valid } =
    if check_valid && not valid then raise Invalidated;
    let ni, nj = Bigarray.Array2.dim2 d - 1, Bigarray.Array2.dim1 d - 1 in
    Format.pp_print_string fmt "[";
    Format.pp_open_vbox fmt 0;
    for i = 0 to ni do

      Format.pp_open_hovbox fmt 4;
      for j = 0 to nj do
        if j > 0 then (
          Format.pp_print_string fmt " ";
          Format.pp_print_cut fmt ();
        );
        Format.fprintf fmt "% -15e" d.{j, i}
      done;
      Format.pp_close_box fmt ();

      if i < ni then (
        Format.pp_print_string fmt ";";
        Format.pp_print_cut fmt ();
      );

    done;
    Format.pp_close_box fmt ();
    Format.pp_print_string fmt "]"

  let ppi ?(start="[") ?(stop="]") ?(sep=";") ?(indent=4) ?(itemsep=" ")
          ?(item=fun f->Format.fprintf f "(%2d,%2d)=% -15e") ()
          fmt { payload = d; valid } =
    if check_valid && not valid then raise Invalidated;
    let ni, nj = Bigarray.Array2.dim2 d - 1, Bigarray.Array2.dim1 d - 1 in
    Format.pp_print_string fmt start;
    Format.pp_open_vbox fmt 0;
    for i = 0 to ni do

      Format.pp_open_hovbox fmt indent;
      for j = 0 to nj do
        if j > 0 then (
          Format.pp_print_string fmt itemsep;
          Format.pp_print_cut fmt ();
        );
        item fmt i j d.{j, i}
      done;
      Format.pp_close_box fmt ();

      if i < ni then (
        Format.pp_print_string fmt sep;
        Format.pp_print_cut fmt ();
      );

    done;
    Format.pp_close_box fmt ();
    Format.pp_print_string fmt stop

  (*
  external c_get : cptr -> int -> int -> float
    = "ml_matrix_dense_get"
  *)

  let get m i j =
    if check_valid && not m.valid then raise Invalidated;
    m.payload.{j, i}

  (* external c_set : cptr -> int -> int -> float -> unit
    = "ml_matrix_dense_set" *)

  let set m i j v =
    if check_valid && not m.valid then raise Invalidated;
    m.payload.{j, i} <- v

  let update m i j f =
    if check_valid && not m.valid then raise Invalidated;
    m.payload.{j, i} <- f m.payload.{j, i}

  external c_scale_add : float -> cptr -> cptr -> unit
    = "sunml_matrix_dense_scale_add"

  let scale_add c m1 m2 =
    if check_valid && not (m1.valid && m2.valid) then raise Invalidated;
    (* compatability checked on C side *)
    c_scale_add c m1.rawptr m2.rawptr

  external c_scale_addi : float -> cptr -> unit
    = "sunml_matrix_dense_scale_addi"

  let scale_addi c m =
    if check_valid && not m.valid then raise Invalidated;
    c_scale_addi c m.rawptr

  external c_matvec
    : cptr -> RealArray.t -> RealArray.t -> unit
    = "sunml_matrix_dense_matvec"

  let matvec m x y =
    if check_valid && not m.valid then raise Invalidated;
    if Sundials_configuration.safe then begin
      let xl = RealArray.length x in
      let yl = RealArray.length y in
      let m, n = Bigarray.Array2.(dim2 m.payload, dim1 m.payload) in
      if n <> xl || m <> yl then raise IncompatibleArguments
    end;
    c_matvec m.rawptr x y

  let set_to_zero { payload = data; valid } =
    if check_valid && not valid then raise Invalidated;
    Bigarray.Array2.fill data 0.0

  let blit ~src ~dst =
    if check_valid && not (src.valid && src.valid) then raise Invalidated;
    Bigarray.Array2.blit src.payload dst.payload

  external c_space : cptr -> int * int
      = "sunml_matrix_dense_space"

  let space { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_space rawptr

  let clone { payload = dataa; valid } =
    if check_valid && not valid then raise Invalidated;
    let m, n = Bigarray.Array2.(dim2 dataa, dim1 dataa) in
    create m n

  let ops = {
    m_clone      = clone;

    m_zero       = set_to_zero;

    m_copy       = (fun src dst -> blit ~src ~dst);

    m_scale_add  = scale_add;

    m_scale_addi = scale_addi;

    m_matvec_setup = None;

    m_matvec     = (fun m x y ->
                      if check_valid && not m.valid then raise Invalidated;
                      c_matvec m.rawptr x y);

    m_space      = space;
  }
end (* }}} *)

module Band = struct (* {{{ *)

  (* Must correspond with sundials_matrix_ml.h:mat_band_dimensions_index *)
  type dimensions = {
      n   : int;
      mu  : int;
      smu : int;
      ml  : int;
    }

  (* Must correspond with sundials_matrix_ml.h:mat_band_data_index *)
  type data = {
    data : RealArray2.data;
    dims : dimensions;
  }
  type cptr

  type t = (data, cptr) matrix_content

  external c_create : dimensions -> t
    = "sunml_matrix_band_create"

  let create ({n; smu; ml; _} as dimensions) =
    if Sundials_configuration.safe then begin
      if n <= 0 then invalid_arg "n";
      if smu < 0 then invalid_arg "smu";
      if ml < 0 then invalid_arg "ml"
    end;
    c_create dimensions

  let unwrap { payload = { data }; _ } = data

  let make dimensions x =
    let { payload = { data }; _ } as r = create dimensions in
    Bigarray.Array2.fill data x;
    r

  let invalidate v =
    if check_valid then v.valid <- false

  let size { payload = { dims = { n } }; valid } =
    if check_valid && not valid then raise Invalidated;
    n, n

  let dims { payload = { dims }; valid; _ } =
    if check_valid && not valid then raise Invalidated;
    dims

  let pp fmt { payload = { data; dims = { n; mu; ml; smu } }; valid } =
    if check_valid && not valid then raise Invalidated;

    Format.pp_print_string fmt "[";
    Format.pp_open_vbox fmt 0;
    for i = 0 to n - 1 do

      Format.pp_open_hovbox fmt 4;
      for j = 0 to n - 1 do
        if j > 0 then (
          Format.pp_print_string fmt " ";
          Format.pp_print_cut fmt ();
        );
        if (i > j + ml) || (j > i + mu)
        then Format.pp_print_string fmt "       ~       "
        else Format.fprintf fmt "% -15e" data.{j, i - j + smu}
      done;
      Format.pp_close_box fmt ();

      if i < n - 1 then (
        Format.pp_print_string fmt ";";
        Format.pp_print_cut fmt ();
      );

    done;
    Format.pp_close_box fmt ();
    Format.pp_print_string fmt "]"

  let ppi ?(start="[") ?(stop="]") ?(sep=";") ?(indent=4) ?(itemsep=" ")
          ?(empty="           ~           ")
          ?(item=fun f->Format.fprintf f "(%2d,%2d)=% -15e") ()
          fmt { payload = { data; dims = { n; mu; ml; smu } }; valid } =
    if check_valid && not valid then raise Invalidated;

    Format.pp_print_string fmt start;
    Format.pp_open_vbox fmt 0;
    for i = 0 to n - 1 do

      Format.pp_open_hovbox fmt indent;
      for j = 0 to n - 1 do
        if j > 0 then (
          Format.pp_print_string fmt itemsep;
          Format.pp_print_cut fmt ();
        );
        if (i > j + ml) || (j > i + mu)
        then Format.pp_print_string fmt empty
        else item fmt i j data.{j, i - j + smu}
      done;
      Format.pp_close_box fmt ();

      if i < n - 1 then (
        Format.pp_print_string fmt sep;
        Format.pp_print_cut fmt ();
      );

    done;
    Format.pp_close_box fmt ();
    Format.pp_print_string fmt stop

  let get m i j =
    let { payload = {data; dims = { smu } }; valid; _ } = m in
    if check_valid && not valid then raise Invalidated;
    data.{j, i - j + smu}

  let set m i j v =
    let { payload={data; dims = { smu } }; valid; _ } = m in
    if check_valid && not valid then raise Invalidated;
    data.{j, i - j + smu} <- v

  let update m i j f =
    let { payload={data; dims = { smu }}; valid; _ } = m in
    if check_valid && not valid then raise Invalidated;
    let k = i - j + smu in
    data.{j, k} <- f data.{j, k}

  external c_scale_add : float -> t -> cptr -> unit
    = "sunml_matrix_band_scale_add"

  let scale_add c m1 m2 =
    if check_valid && not (m1.valid && m2.valid) then raise Invalidated;
    c_scale_add c m1 m2.rawptr

  external c_scale_addi : float -> cptr -> unit
    = "sunml_matrix_band_scale_addi"

  let scale_addi c { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_scale_addi c rawptr

  external c_matvec
    : cptr -> RealArray.t -> RealArray.t -> unit
    = "sunml_matrix_band_matvec"

  let matvec m x y =
    if check_valid && not m.valid then raise Invalidated;
    if Sundials_configuration.safe then begin
      let xl = RealArray.length x in
      let yl = RealArray.length y in
      let { dims = { n }; _ } = m.payload in
      if n <> yl || n <> xl then raise IncompatibleArguments
    end;
    c_matvec m.rawptr x y

  let set_to_zero { payload = { data }; valid; _ } =
    if check_valid && not valid then raise Invalidated;
    Bigarray.Array2.fill data 0.0

  external c_copy : cptr -> t -> unit
    = "sunml_matrix_band_copy"

  let blit ~src ~dst =
    if check_valid && not (src.valid && dst.valid) then raise Invalidated;
    c_copy src.rawptr dst

  let space { payload = { dims = { n; smu; ml; _ } }; valid } =
    if check_valid && not valid then raise Invalidated;
    (n * (smu + ml + 1), 7 + n)

  let clone { payload = { dims; _ }; valid } =
    if check_valid && not valid then raise Invalidated;
    create dims

  let ops = {
    m_clone      = clone;

    m_zero       = set_to_zero;

    m_copy       = (fun src dst -> blit ~src ~dst);

    m_scale_add  = scale_add;

    m_scale_addi = scale_addi;

    m_matvec_setup = None;

    m_matvec     = (fun m x y ->
                      if check_valid && not m.valid then raise Invalidated;
                      c_matvec m.rawptr x y);

    m_space      = space;
  }
end (* }}} *)

module Sparse = struct (* {{{ *)
  type csc
  type csr

  (* See MAT_TO/FROM_SFORMAT(x) in sundials_matrix_ml.h *)
  type _ sformat =
    | CSC : csc sformat (* Must correpond with Sundials CSC_MAT = 0 constant. *)
    | CSR : csr sformat (* Must correpond with Sundials CSR_MAT = 1 constant. *)

  (* The payload field cannot be relied upon for Sundials < 3.0.0 *)
  let unsafe_content =
    match Config.sundials_version with
    | 2,_,_ -> true
    | _ -> false

  type index_array =
    (Index.t, Index.index_elt, Bigarray.c_layout) Bigarray.Array1.t

  (* Must correspond with sundials_matrix_ml.h:mat_sparse_data_index *)
  type 's data = {
    idxvals : index_array;
    idxptrs : index_array;
    data    : RealArray.t;
    sformat : 's sformat;
  }

  type cptr

  type 's t = ('s data, cptr) matrix_content

  external c_create : int -> int -> int -> 's sformat -> 's t
    = "sunml_matrix_sparse_create"

  let make (type s) (sformat : s sformat) m n nnz =
    if Sundials_configuration.safe then begin
      if m <= 0 then invalid_arg "number of rows";
      if n <= 0 then invalid_arg "number of columns";
      if nnz < 0 then invalid_arg "number of non-zeros"
    end;
    if Sundials_configuration.safe then
      (match Config.sundials_version, sformat with
       | (2,v,_), _ when v < 6 ->
           raise Config.NotImplementedBySundialsVersion
       | (2,v,_), CSR when v < 7 ->
           raise Config.NotImplementedBySundialsVersion
       | _ -> ());
    c_create m n nnz sformat

  external c_from_dense : 's sformat -> Dense.cptr -> float -> 's t
    = "sunml_matrix_sparse_from_dense"

  let from_dense (type s) (sformat : s sformat)
        droptol { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    if Sundials_configuration.safe && droptol < 0.0 then invalid_arg "droptol";
    if Sundials_configuration.safe then
      (match Config.sundials_version, sformat with
       | (2,v,_), _ when v < 6 ->
           raise Config.NotImplementedBySundialsVersion
       | (2,v,_), CSR when v < 7 ->
           raise Config.NotImplementedBySundialsVersion
       | _ -> ());
    c_from_dense sformat rawptr droptol

  external c_from_band : 's sformat -> Band.cptr -> float -> 's t
    = "sunml_matrix_sparse_from_band"

  let from_band (type s) (sformat : s sformat)
        droptol { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    if Sundials_configuration.safe && droptol < 0.0 then invalid_arg "droptol";
    if Sundials_configuration.safe then
      (match Config.sundials_version, sformat with
       | (2,v,_), _ when v < 6 ->
           raise Config.NotImplementedBySundialsVersion
       | (2,v,_), CSR when v < 7 ->
           raise Config.NotImplementedBySundialsVersion
       | _ -> ());
    c_from_band sformat rawptr droptol

  let sformat { payload = { sformat } } = sformat

  let is_csc (type s) (mat : s t) =
    match sformat mat with
    | CSC -> true
    | CSR -> false

  external c_rewrap : 's t -> 's data
    = "sunml_matrix_sparse_rewrap"

  let unwrap ({ payload; _ } as m) =
    let { idxvals; idxptrs; data } =
      if unsafe_content then c_rewrap m else payload
    in
    idxvals, idxptrs, data

  external c_resize : 's t -> int -> bool -> unit
    = "sunml_matrix_sparse_resize"

  let resize ?nnz a =
    let nnz = match nnz with Some x -> x | None -> 0 in
    c_resize a nnz true

  let invalidate v =
    if check_valid then v.valid <- false

  external c_size : cptr -> int * int
    = "sunml_matrix_sparse_size"

  let size { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_size rawptr

  external c_dims : cptr -> int * int
    = "sunml_matrix_sparse_dims"

  let dims { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_dims rawptr

  external c_set_idx : cptr -> int -> int -> unit
    = "sunml_matrix_sparse_set_idx"

  external c_get_idx : cptr -> int -> int
    = "sunml_matrix_sparse_get_idx"

  let set_col m j idx =
    let { payload = { idxptrs }; rawptr; valid } = m in
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_set_idx rawptr j idx
    else idxptrs.{j} <- Index.of_int idx

  let get_col m j =
    let { payload = { idxptrs }; rawptr; valid } = m in
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_get_idx rawptr j
    else Index.to_int idxptrs.{j}

  let set_row = set_col
  let get_row = get_col

  external c_set_data : cptr -> int -> float -> unit
    = "sunml_matrix_sparse_set_data"

  external c_get_data : cptr -> int -> float
    = "sunml_matrix_sparse_get_data"

  external c_set_val : cptr -> int -> int -> unit
    = "sunml_matrix_sparse_set_val"

  external c_get_val : cptr -> int -> int
    = "sunml_matrix_sparse_get_val"

  let set m idx i v =
    let{ payload = { idxvals; data }; rawptr; valid } = m in
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then (c_set_val rawptr idx i; c_set_data rawptr idx v)
    else (idxvals.{idx} <- Index.of_int i; data.{idx} <- v)

  let set_data m idx v =
    let { payload = { data }; rawptr; valid } = m in
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_set_data rawptr idx v
    else data.{idx} <- v

  let set_rowval m idx i =
    let { payload = { idxvals }; rawptr; valid } = m in
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_set_val rawptr idx i
    else idxvals.{idx} <- Index.of_int i

  let set_colval = set_rowval

  let get m idx =
    let { payload = { idxvals; data }; rawptr; valid } = m in
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_get_val rawptr idx, c_get_data rawptr idx
    else Index.to_int idxvals.{idx}, data.{idx}

  let get_rowval m idx =
    let { payload = { idxvals }; rawptr; valid } = m in
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_get_val rawptr idx
    else Index.to_int idxvals.{idx}

  let get_colval = get_rowval

  let get_data m idx =
    let { payload = { data }; rawptr; valid } = m in
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_get_data rawptr idx
    else data.{idx}

  let pp (type s) fmt (mat : s t) =
    if check_valid && not mat.valid then raise Invalidated;
    let _, n = size mat in
    let end_row = ref false in
    let idxrole = match mat.payload.sformat with
                  | CSC -> "col"
                  | CSR -> "row" in

    Format.pp_print_string fmt "[";
    Format.pp_open_vbox fmt 0;
    for j = 0 to n - 1 do
      let p, np = get_col mat j, get_col mat (j + 1) in
      if p < np then begin
        if !end_row then begin
          Format.pp_print_string fmt ";";
          Format.pp_close_box fmt ();
          Format.pp_print_cut fmt ()
        end;
        Format.pp_open_hovbox fmt 0;
        Format.fprintf fmt "%s %2d: " idxrole j;
        Format.pp_print_space fmt ();

        for i = p to np - 1 do
          if i > p then Format.pp_print_space fmt ();
          let r, v = get mat i in
          Format.fprintf fmt "%2d=% -15e" r v
        done;
        end_row := true
      end;
    done;
    if !end_row then Format.pp_close_box fmt ();
    Format.pp_close_box fmt ();
    Format.pp_print_string fmt "]"

  let ppi ?(start="[") ?(stop="]") ?(sep=";") ?(indent=0) ?(itemsep=" ")
          ?(rowcol=fun f->Format.fprintf f "%2d: ")
          ?(item=fun f->Format.fprintf f "%2d=% -15e") ()
          fmt mat =
    if check_valid && not mat.valid then raise Invalidated;
    let _, n = size mat in
    let end_row = ref false in

    Format.pp_print_string fmt start;
    Format.pp_open_vbox fmt 0;
    for j = 0 to n - 1 do
      let p, np = get_col mat j, get_col mat (j + 1) in
      if p < np then begin
        if !end_row then begin
          Format.pp_print_string fmt sep;
          Format.pp_close_box fmt ();
          Format.pp_print_cut fmt ()
        end;
        Format.pp_open_hovbox fmt indent;
        rowcol fmt j;
        Format.pp_print_space fmt ();

        for i = p to np - 1 do
          if i > p then (
            Format.pp_print_string fmt itemsep;
            Format.pp_print_cut fmt ();
          );
          let r, v = get mat i in
          item fmt r v
        done;
        end_row := true
      end;
    done;
    if !end_row then Format.pp_close_box fmt ();
    Format.pp_close_box fmt ();
    Format.pp_print_string fmt stop

  external c_scale_add : float -> 's t -> cptr -> unit
    = "sunml_matrix_sparse_scale_add"

  let scale_add c m1 m2 =
    if check_valid && not (m1.valid && m2.valid) then raise Invalidated;
    if Sundials_configuration.safe && c_size m1.rawptr <> c_size m2.rawptr
    then raise IncompatibleArguments;
    c_scale_add c m1 m2.rawptr

  external c_scale_addi : float -> 's t -> unit
    = "sunml_matrix_sparse_scale_addi"

  let scale_addi c ({ valid } as a) =
    if check_valid && not valid then raise Invalidated;
    c_scale_addi c a

  external c_matvec
    : cptr -> RealArray.t -> RealArray.t -> unit
    = "sunml_matrix_sparse_matvec"

  let matvec m x y =
    if check_valid && not m.valid then raise Invalidated;
    if Sundials_configuration.safe then begin
      let xl = RealArray.length x in
      let yl = RealArray.length y in
      let m, n = c_size m.rawptr in
      if n <> xl || m <> yl then raise IncompatibleArguments
    end;
    c_matvec m.rawptr x y

  external c_set_to_zero : cptr -> unit
    = "sunml_matrix_sparse_set_to_zero"

  let set_to_zero { payload = { idxvals; idxptrs; data }; rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_set_to_zero rawptr
    else (Bigarray.Array1.fill idxvals Index.zero;
          Bigarray.Array1.fill idxptrs Index.zero;
          Bigarray.Array1.fill data 0.0)

  external c_copy : cptr -> 's t -> unit
      = "sunml_matrix_sparse_copy"

  let blit ~src ~dst =
    if check_valid && not (src.valid && dst.valid) then raise Invalidated;
    if Sundials_configuration.safe && c_size src.rawptr <> c_size dst.rawptr
    then raise IncompatibleArguments;
    c_copy src.rawptr dst

  external tocsr : cptr -> csr t = "sunml_matrix_sparse_tocsr"
  external tocsc : cptr -> csc t = "sunml_matrix_sparse_tocsc"

  let copy_to_csr { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    tocsr rawptr

  let copy_to_csc { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    tocsc rawptr

  external c_space : cptr -> int * int
      = "sunml_matrix_sparse_space"

  let space { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_space rawptr

  let clone { rawptr; valid; payload = { sformat = fmta; _ } } =
    if check_valid && not valid then raise Invalidated;
    let m, n = c_size rawptr in
    let nnz, _ = c_dims rawptr in
    c_create m n nnz fmta

  let ops = {
    m_clone      = clone;

    m_zero       = set_to_zero;

    m_copy       = (fun src dst -> blit ~src ~dst);

    m_scale_add  = scale_add;

    m_scale_addi = scale_addi;

    m_matvec_setup = None;

    m_matvec     = (fun m x y ->
                      if check_valid && not m.valid then raise Invalidated;
                      c_matvec m.rawptr x y);

    m_space      = space;
  }
end (* }}} *)

type lint_array = LintArray.t
type real_array = RealArray.t

module ArrayDense = struct (* {{{ *)
  type t = RealArray2.t

  let make = RealArray2.make
  let create = RealArray2.create
  let size = RealArray2.size
  let pp = RealArray2.pp
  let unwrap = RealArray2.unwrap
  let get = RealArray2.get
  let set = RealArray2.set

  let update a i j f = set a i j (f (get a i j))

  let set_to_zero x = Bigarray.Array2.fill (RealArray2.unwrap x) 0.0

  let blit = RealArray2.blit

  external scale : float -> t -> unit
      = "sunml_arraydensematrix_scale"

  external add_identity : t -> unit
      = "sunml_arraydensematrix_add_identity"

  external matvec' : t -> real_array -> real_array -> unit
      = "sunml_arraydensematrix_matvec"

  let matvec a x y = matvec' a x y

  external getrf : t -> lint_array -> unit
      = "sunml_arraydensematrix_getrf"

  external getrs : t -> lint_array -> real_array -> unit
      = "sunml_arraydensematrix_getrs"

  external getrs' : t -> lint_array -> real_array -> int -> unit
      = "sunml_arraydensematrix_getrs_off"

  external potrf : t -> unit
      = "sunml_arraydensematrix_potrf"

  external potrs : t -> real_array -> unit
      = "sunml_arraydensematrix_potrs"

  external geqrf : t -> real_array -> real_array -> unit
      = "sunml_arraydensematrix_geqrf"

  external ormqr'
      : t -> (real_array * real_array * real_array * real_array) -> unit
      = "sunml_arraydensematrix_ormqr"

  let ormqr ~a ~beta ~v ~w ~work = ormqr' a (beta, v, w, work)

  let clone a =
    let m, n = RealArray2.size a in
    RealArray2.create m n

  let scale_addi c a =
    let ad = RealArray2.unwrap a in
    let an, am = Bigarray.Array2.(dim1 ad, dim2 ad) in
    for j = 0 to am - 1 do
      for i = 0 to an - 1 do
        ad.{i, j} <- c *. ad.{i, j} +. (if i = j then 1. else 0.)
      done
    done

  let scale_add c a b =
    let ad, bd = RealArray2.unwrap a, RealArray2.unwrap b in
    let an, am = Bigarray.Array2.(dim1 ad, dim2 ad) in
    let bn, bm = Bigarray.Array2.(dim1 bd, dim2 bd) in
    if Sundials_configuration.safe && (am <> bm || an <> bn)
    then raise IncompatibleArguments;
    for i = 0 to an - 1 do
      for j =0 to am - 1 do
        ad.{i, j} <- c *. ad.{i, j} +. bd.{i, j}
      done
    done

  let space a =
    let m, n = RealArray2.size a in
    (m * n, 3)

  let ops = {
    m_clone      = clone;

    m_zero       = set_to_zero;

    m_copy       = (fun src dst -> blit ~src ~dst);

    m_scale_add  = scale_add;

    m_scale_addi = scale_addi;

    m_matvec_setup = None;

    m_matvec     = matvec;

    m_space      = space;
  }
end (* }}} *)

module ArrayBand = struct (* {{{ *)

  type smu = int
  type mu = int
  type ml = int

  type t = RealArray2.t * (smu * mu * ml)

  let make ((smu, _, ml) as dims) n v =
    (RealArray2.make (smu + ml + 1) n v, dims)

  let create ((smu, _, ml) as dims) n =
    (RealArray2.create (smu + ml + 1) n, dims)

  let size (a, _) = let (_, m) = RealArray2.size a in m, m
  let dims (_, dims) = dims

  let get (a, (smu, _, _)) i j =
    RealArray2.get a (i - j + smu) j

  let set (a, (smu, _, _)) i j v =
    RealArray2.set a (i - j + smu) j v

  let update (a, (smu, _, _)) i j f =
    let k = i - j + smu in
    RealArray2.set a k j (f (RealArray2.get a k j))

  let pp fmt (a, (smu, mu, ml)) =
    let data = RealArray2.unwrap a in
    let _, n = RealArray2.size a in
    Format.pp_print_string fmt "[";
    Format.pp_open_vbox fmt 0;
    for i = 0 to n - 1 do

      Format.pp_open_hovbox fmt 4;
      for j = 0 to n - 1 do
        if j > 0 then (
          Format.pp_print_string fmt " ";
          Format.pp_print_cut fmt ();
        );
        if (i > j + ml) || (j > i + mu)
        then Format.pp_print_string fmt "       ~       "
        else Format.fprintf fmt "% -15e" data.{j, i - j + smu}
      done;
      Format.pp_close_box fmt ();

      if i < n - 1 then (
        Format.pp_print_string fmt ";";
        Format.pp_print_cut fmt ();
      );

    done;
    Format.pp_close_box fmt ();
    Format.pp_print_string fmt "]"

  let ppi ?(start="[") ?(stop="]") ?(sep=";") ?(indent=4) ?(itemsep=" ")
          ?(empty="           ~           ")
          ?(item=fun f->Format.fprintf f "(%2d,%2d)=% -15e") ()
          fmt (a, (smu, mu, ml)) =
    let data = RealArray2.unwrap a in
    let _, n = RealArray2.size a in

    Format.pp_print_string fmt start;
    Format.pp_open_vbox fmt 0;
    for i = 0 to n - 1 do

      Format.pp_open_hovbox fmt indent;
      for j = 0 to n - 1 do
        if j > 0 then (
          Format.pp_print_string fmt itemsep;
          Format.pp_print_cut fmt ();
        );
        if (i > j + ml) || (j > i + mu)
        then Format.pp_print_string fmt empty
        else item fmt i j data.{j, i - j + smu}
      done;
      Format.pp_close_box fmt ();

      if i < n - 1 then (
        Format.pp_print_string fmt sep;
        Format.pp_print_cut fmt ();
      );

    done;
    Format.pp_close_box fmt ();
    Format.pp_print_string fmt stop

  let unwrap (a, _) = RealArray2.unwrap a

  external copy'
    : RealArray2.t -> RealArray2.t
      -> int * int * int * int -> unit
    = "sunml_arraybandmatrix_copy"

  let set_to_zero (x, _) =
    Bigarray.Array2.fill (RealArray2.unwrap x) 0.0

  let blit ~src:(a, (a_smu, a_mu, a_ml)) ~dst:(b, (b_smu, _, _))
      = copy' a b (a_smu, b_smu, a_mu, a_ml)

  external scale' : float -> RealArray2.t -> int * int * int -> unit
      = "sunml_arraybandmatrix_scale"

  let scale c (a, dims) = scale' c a dims

  external c_add_identity : smu -> RealArray2.t -> unit
      = "sunml_arraybandmatrix_add_identity"

  let add_identity (x, (smu, _, _)) = c_add_identity smu x

  external c_matvec
    : RealArray2.t -> int * int * int
      -> real_array -> real_array -> unit
    = "sunml_arraybandmatrix_matvec"

  let matvec (a, dims) x y = c_matvec a dims x y

  external gbtrf'
    : RealArray2.t -> int * int * int -> lint_array -> unit
    = "sunml_arraybandmatrix_gbtrf"

  let gbtrf (a, dims) p = gbtrf' a dims p

  external gbtrs'
    : RealArray2.t -> int * int * int -> lint_array -> real_array -> unit
    = "sunml_arraybandmatrix_gbtrs"

  let gbtrs (a, dims) p b = gbtrs' a dims p b

  let clone (a, dims) =
    let m, n = RealArray2.size a in
    (RealArray2.make m n 0.0, dims)

  let scale_addi c (a, (smu, mu, ml)) =
    let ad = RealArray2.unwrap a in
    let _, am = RealArray2.size a in
    for j = 0 to am - 1 do
      for i = smu - mu to smu + ml do
        ad.{j, i} <- c *. ad.{j, i}
      done;
      ad.{j, smu} <- ad.{j, smu} +. 1.
    done

  let scale_add c (a, (smu_a, mu_a, ml_a)) (b, (smu_b, mu_b, ml_b)) =
    let ad, bd = RealArray2.unwrap a, RealArray2.unwrap b in
    let (_, am), (_, bm) = RealArray2.size a, RealArray2.size b in
    if Sundials_configuration.safe && (am <> bm || mu_b > mu_a || ml_b > ml_a)
    then raise IncompatibleArguments;
    for j = 0 to bm - 1 do
      for i = - mu_b to ml_b do
        ad.{j, i + smu_a} <- c *. ad.{j, i + smu_a} +. bd.{j, i + smu_b}
      done;
    done

  let space (a, (smu, _, ml)) =
    let _, n = RealArray2.size a in
    (n * (smu + ml + 1), 7 + n)
    (* 3 integer variables + 2 integer dimensions *)

  let ops = {
    m_clone      = clone;

    m_zero       = set_to_zero;

    m_copy       = (fun src dst -> blit ~src ~dst);

    m_scale_add  = scale_add;

    m_scale_addi = scale_addi;

    m_matvec_setup = None;

    m_matvec     = matvec;

    m_space      = space;
  }

end (* }}} *)

type standard
type custom

(* Must correspond with sundials_matrix_ml.h:mat_matrix_id_tag *)
type (_,_,_,_) id =
  | Dense : (standard, Dense.t, Nvector_serial.data, [>Nvector_serial.kind]) id
  | Band  : (standard, Band.t, Nvector_serial.data, [>Nvector_serial.kind]) id
  | SparseCSC : (standard, Sparse.csc Sparse.t, Nvector_serial.data, [>Nvector_serial.kind]) id
  | SparseCSR : (standard, Sparse.csr Sparse.t, Nvector_serial.data, [>Nvector_serial.kind]) id
  | Custom : (custom, 'm, 'nd, 'nk) id
  | ArrayDense : (custom, ArrayDense.t, RealArray.t, 'nk) id
  | ArrayBand  : (custom, ArrayBand.t, RealArray.t, 'nk) id

type cmat

(* Must correspond with sundials_matrix_ml.h:mat_matrix_index *)
type [@warning "-69"] ('k, 'm, 'nd, 'nk) t = {
  payload : 'm;
  rawptr  : cmat;
  id      : ('k, 'm, 'nd, 'nk) id;
  mat_ops : ('m, 'nd) matrix_ops;
}

(* The labels must also correspond with sundials_matrix_ml.h:mat_matrix_id_tag *)
type (_, _) any =
  | ADense : (standard, Dense.t, Nvector_serial.data, [>Nvector_serial.kind]) t
             -> (Nvector_serial.data, [>Nvector_serial.kind]) any
  | ABand  : (standard, Band.t, Nvector_serial.data, [>Nvector_serial.kind]) t
             -> (Nvector_serial.data, [>Nvector_serial.kind]) any
  | ASparseCSC : (standard, Sparse.csc Sparse.t, Nvector_serial.data, [>Nvector_serial.kind]) t
                 -> (Nvector_serial.data, [>Nvector_serial.kind]) any
  | ASparseCSR : (standard, Sparse.csr Sparse.t, Nvector_serial.data, [>Nvector_serial.kind]) t
                 -> (Nvector_serial.data, [>Nvector_serial.kind]) any
  | ACustom : (custom, 'm, 'nd, 'nk) t -> ('nd, 'nk) any
  | AArrayDense : (custom, ArrayDense.t, RealArray.t, 'nk) t -> (RealArray.t, 'nk) any
  | AArrayBand  : (custom, ArrayBand.t, RealArray.t, 'nk) t -> (RealArray.t, 'nk) any

type 'nk dense =
  (standard, Dense.t, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

type 'nk band =
  (standard, Band.t, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

type ('s, 'nk) sparse =
  (standard, 's Sparse.t, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

type 'nk arraydense = (custom, ArrayDense.t, RealArray.t, 'nk) t

type 'nk arrayband = (custom, ArrayBand.t, RealArray.t, 'nk) t

external c_wrap :
     ('k, 'm, 'nd, 'nk) id
  -> 'content_cptr
  -> 'm
  -> bool
  -> Context.t
  -> cmat
  = "sunml_matrix_wrap"

let wrap_dense ?context (data : Dense.t) =
  let ctx = Sundials_impl.Context.get context in
  {
    payload = data;
    rawptr  = c_wrap Dense data.rawptr data false ctx;
    id      = Dense;
    mat_ops = Dense.ops;
  }

let dense ?context ?m ?(i=0.0) n =
  let m = match m with Some m -> m | None -> n in
  wrap_dense ?context (Dense.make m n i)

let wrap_band ?context (data : Band.t) =
  let ctx = Sundials_impl.Context.get context in
  {
    payload = data;
    rawptr  = c_wrap Band data.rawptr data false ctx;
    id      = Band;
    mat_ops = Band.ops;
  }

let band ?context ?(mu=2) ?smu ?ml ?(i=0.0) n =
  let ml  = match ml with Some ml -> ml | None -> mu in
  let smu = match smu with Some smu -> smu | None -> mu+ml in
  wrap_band ?context (Band.(make { n; mu; smu; ml } i))

let wrap_sparse (type sf) ?context (data : sf Sparse.t) =
  (match Config.sundials_version with
   | (2,v,_) when v < 6 -> raise Config.NotImplementedBySundialsVersion
   | _ -> ());
  let sparse_id : (standard, sf Sparse.t, Nvector_serial.data, [>Nvector_serial.kind]) id =
    match Sparse.sformat data with
    | Sparse.CSC -> SparseCSC
    | Sparse.CSR -> SparseCSR
  in
  let ctx = Sundials_impl.Context.get context in
  {
    payload = data;
    rawptr  = c_wrap sparse_id data.rawptr data false ctx;
    id      = sparse_id;
    mat_ops = Sparse.ops;
  }

let sparse_csc ?context ?m ?nnz n =
  let m = match m with Some m -> m | None -> n in
  let nnz = match nnz with Some nnz -> nnz | None -> n / 10 in
  wrap_sparse ?context Sparse.(make CSC m n nnz)

let sparse_csr ?context ?m ?nnz n =
  let m = match m with Some m -> m | None -> n in
  let nnz = match nnz with Some nnz -> nnz | None -> n / 10 in
  wrap_sparse ?context Sparse.(make CSR m n nnz)

let wrap_custom ops ?context data =
  let ctx = Sundials_impl.Context.get context in
  {
    payload = data;
    rawptr  = c_wrap Custom (ops, Custom) data (ops.m_matvec_setup <> None) ctx;
    id      = Custom;
    mat_ops = ops;
  }

let wrap_arraydense ?context data =
  let ctx = Sundials_impl.Context.get context in
  {
    payload = data;
    rawptr  = c_wrap ArrayDense (ArrayDense.ops, ArrayDense) data false ctx;
    id      = ArrayDense;
    mat_ops = ArrayDense.ops;
  }

let arraydense ?context ?m ?(i=0.0) n =
  let m = match m with Some m -> m | None -> n in
  wrap_arraydense ?context (ArrayDense.make m n i)

let wrap_arrayband ?context data =
  let ctx = Sundials_impl.Context.get context in
  {
    payload = data;
    rawptr  = c_wrap ArrayBand (ArrayBand.ops, ArrayBand) data false ctx;
    id      = ArrayBand;
    mat_ops = ArrayBand.ops;
  }

let arrayband ?context ?mu ?smu ?ml ?(i=0.0) n =
  let mu = match mu, smu with
           | None, None  -> 2
           | Some mu, _  -> mu
           | _, Some smu -> smu
  in
  let ml  = match ml  with Some ml  -> ml  | None -> mu in
  let smu = match smu with Some smu -> smu | None -> mu+ml in
  wrap_arrayband ?context (ArrayBand.make (smu, mu, ml) n i)

let get_ops { mat_ops } = mat_ops

let get_id { id } = id

let unwrap { payload } = payload

external c_scale_add
  : float -> ('k, 'm, 'nd, 'nk) t -> ('k, 'm, 'nd, 'nk) t -> unit
  = "sunml_matrix_scale_add"

let scale_add c ({ payload = a; mat_ops = { m_scale_add } } as m1)
                ({ payload = b } as m2) =
  match Config.sundials_version with
  | 2,_,_ -> m_scale_add c a b
  | _ -> c_scale_add c m1 m2

external c_scale_addi : float -> ('k, 'm, 'nd, 'nk) t -> unit
  = "sunml_matrix_scale_addi"

let scale_addi c ({ payload = a; mat_ops = { m_scale_addi } } as m) =
  match Config.sundials_version with
  | 2,_,_ -> m_scale_addi c a
  | _ -> c_scale_addi c m

external c_matvecsetup : ('k, 'm, 'nd, 'nk) t -> unit
  = "sunml_matrix_matvecsetup"

let matvec_setup m =
  if Sundials_impl.Version.lt500
    then raise Config.NotImplementedBySundialsVersion;
  c_matvecsetup m

external c_matvec
  : ('k, 'm, 'nd, 'nk) t -> ('nd, 'nk) Nvector.t -> ('nd, 'nk) Nvector.t -> unit
  = "sunml_matrix_matvec"

let matvec ({ payload = a; mat_ops = { m_matvec } } as m) x y =
  match Config.sundials_version with
  | 2,_,_ -> m_matvec a (Nvector.unwrap x) (Nvector.unwrap y)
  | _ -> c_matvec m x y

external c_zero : ('k, 'm, 'nd, 'nk) t -> unit
  = "sunml_matrix_zero"

let set_to_zero ({ payload = a; mat_ops = { m_zero } } as m) =
  match Config.sundials_version with
  | 2,_,_ -> m_zero a
  | _ -> c_zero m

external c_copy : ('k, 'm, 'nd, 'nk) t -> ('k, 'm, 'nd, 'nk) t -> unit
  = "sunml_matrix_copy"

let blit ~src:({ payload = a; mat_ops = { m_copy } } as m1)
         ~dst:({ payload = b } as m2) =
  match Config.sundials_version with
  | 2,_,_ -> m_copy a b
  | _ -> c_copy m1 m2

external c_space : ('k, 'm, 'nd, 'nk) t -> int * int
    = "sunml_matrix_space"

let space ({ payload = a; mat_ops = { m_space } } as m) =
  match Config.sundials_version with
  | 2,_,_ -> m_space a
  | _ -> c_space m

external print_dense : 'nk dense -> Logfile.t -> unit
    = "sunml_matrix_print_dense"

external print_band : 'nk band -> Logfile.t -> unit
    = "sunml_matrix_print_band"

external print_sparse : ('s, 'nk) sparse -> Logfile.t -> unit
    = "sunml_matrix_print_sparse"

let pp (type k m nd nk) fmt ({ id; payload } : (k, m, nd, nk) t) =
  match id with
  | Dense  -> Dense.pp fmt payload
  | Band   -> Band.pp fmt payload
  | SparseCSC -> Sparse.pp fmt payload
  | SparseCSR -> Sparse.pp fmt payload
  | Custom -> Format.pp_print_string fmt "<custom matrix>"
  | ArrayDense -> ArrayDense.pp fmt payload
  | ArrayBand  -> ArrayBand.pp fmt payload

let pp_any (type nd nk) fmt (am : (nd, nk) any) =
  match am with
  | ADense m -> pp fmt m
  | ABand m -> pp fmt m
  | ASparseCSC m -> pp fmt m
  | ASparseCSR m -> pp fmt m
  | ACustom m -> pp fmt m
  | AArrayDense m -> pp fmt m
  | AArrayBand  m -> pp fmt m

(* Let C code know about some of the values in this module.  *)
external c_init_module
  : exn array
    -> (  (Dense.t, Sundials.RealArray.t) matrix_ops
        * (Band.t, Sundials.RealArray.t) matrix_ops
        * ('a Sparse.t, Sundials.RealArray.t) matrix_ops)
    -> unit =
  "sunml_mat_init_module"

let _ =
  c_init_module
    (* Exceptions must be listed in the same order as
       mat_exn_index.  *)
    [|Invalidated;
      IncompatibleArguments;
      ZeroDiagonalElement 0;
    |]
    (Dense.ops, Band.ops, Sparse.ops)

