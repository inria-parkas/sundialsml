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
type ('data, 'cptr) matrix_content = {
  mutable payload : 'data;
  rawptr  : 'cptr;
  mutable valid : bool;
}

type 'k serial_nvector
  = (Nvector_serial.data, [>Nvector_serial.kind] as 'k) Nvector.t

(* Must correspond with matrix_ml.h:mat_matrix_ops_index *)
type ('m, 'd) matrix_ops = {
  m_clone      : 'm -> 'm;

  m_zero       : 'm -> unit;

  m_copy       : 'm -> 'm -> unit;

  m_scale_add  : float -> 'm -> 'm -> unit;

  m_scale_addi : float -> 'm -> unit;

  m_matvec     : 'm -> 'd -> 'd -> unit;

  m_space      : 'm -> int * int;
}

module Dense = struct (* {{{ *)

  type data = RealArray2.data
  type cptr

  type t = (data, cptr) matrix_content

  external c_create : int -> int -> t
    = "sunml_matrix_dense_create"

  let create i j =
    if Sundials_configuration.safe then begin
      if i <= 0 then invalid_arg "i";
      if j <= 0 then invalid_arg "j"
    end;
    c_create i j

  let unwrap { payload } = payload

  let make m n v =
    let { payload } as r = create m n in
    Bigarray.Array2.fill payload 0.0;
    r

  let invalidate v =
    if check_valid then v.valid <- false

  let size { payload; rawptr; valid } =
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

  let get { payload; valid } i j =
    if check_valid && not valid then raise Invalidated;
    payload.{j, i}

  (* external c_set : cptr -> int -> int -> float -> unit
    = "ml_matrix_dense_set" *)

  let set { payload; valid } i j v =
    if check_valid && not valid then raise Invalidated;
    payload.{j, i} <- v

  let update { payload; valid } i j f =
    if check_valid && not valid then raise Invalidated;
    payload.{j, i} <- f payload.{j, i}

  external c_scale_add : float -> cptr -> cptr -> unit
    = "sunml_matrix_dense_scale_add"

  let scale_add c { rawptr = ptr1; valid = valid1 }
                  { rawptr = ptr2; valid = valid2 } =
    if check_valid && not (valid1 && valid2) then raise Invalidated;
    (* compatability checked on C side *)
    c_scale_add c ptr1 ptr2

  external c_scale_addi : float -> cptr -> unit
    = "sunml_matrix_dense_scale_addi"

  let scale_addi c { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_scale_addi c rawptr

  external c_matvec
    : cptr -> RealArray.t -> RealArray.t -> unit
    = "sunml_matrix_dense_matvec"

  let matvec { rawptr; payload; valid } x y =
    if check_valid && not valid then raise Invalidated;
    if Sundials_configuration.safe then begin
      let xl = RealArray.length x in
      let yl = RealArray.length y in
      let m, n = Bigarray.Array2.(dim2 payload, dim1 payload) in
      if n <> xl || m <> yl then raise IncompatibleArguments
    end;
    c_matvec rawptr x y

  let set_to_zero { payload = data; rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    Bigarray.Array2.fill data 0.0

  let blit { payload = data1; valid = valid1 }
           { payload = data2; valid = valid2 } =
    if check_valid && not (valid1 && valid2) then raise Invalidated;
    Bigarray.Array2.blit data1 data2

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

    m_copy       = blit;

    m_scale_add  = scale_add;

    m_scale_addi = scale_addi;

    m_matvec     = (fun { rawptr; valid } x y ->
                      if check_valid && not valid then raise Invalidated;
                      c_matvec rawptr x y);

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

  let create ({n; mu; smu; ml} as dimensions) =
    if Sundials_configuration.safe then begin
      if n <= 0 then invalid_arg "n";
      if smu < 0 then invalid_arg "smu";
      if ml < 0 then invalid_arg "ml"
    end;
    c_create dimensions

  let unwrap { payload = { data }; rawptr } = data

  let make dimensions x =
    let { payload = { data }; rawptr } as r = create dimensions in
    Bigarray.Array2.fill data x;
    r

  let invalidate v =
    if check_valid then v.valid <- false

  let size { payload = { dims = { n } }; valid } =
    if check_valid && not valid then raise Invalidated;
    n, n

  let dims { payload = { dims }; rawptr; valid } =
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

  let get { payload = {data; dims = { smu } }; rawptr; valid } i j =
    if check_valid && not valid then raise Invalidated;
    data.{j, i - j + smu}

  let set { payload={data; dims = { smu } }; rawptr; valid } i j v =
    if check_valid && not valid then raise Invalidated;
    data.{j, i - j + smu} <- v

  let update { payload={data; dims = { smu }}; rawptr; valid } i j f =
    if check_valid && not valid then raise Invalidated;
    let k = i - j + smu in
    data.{j, k} <- f data.{j, k}

  external c_scale_add : float -> t -> cptr -> unit
    = "sunml_matrix_band_scale_add"

  let scale_add c ({ valid = valid1 } as m1)
                  { rawptr = ptr2; valid = valid2 } =
    if check_valid && not (valid1 && valid2) then raise Invalidated;
    c_scale_add c m1 ptr2

  external c_scale_addi : float -> cptr -> unit
    = "sunml_matrix_band_scale_addi"

  let scale_addi c { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_scale_addi c rawptr

  external c_matvec
    : cptr -> RealArray.t -> RealArray.t -> unit
    = "sunml_matrix_band_matvec"

  let matvec { rawptr; payload = { data; dims = { n } }; valid } x y =
    if check_valid && not valid then raise Invalidated;
    if Sundials_configuration.safe then begin
      let xl = RealArray.length x in
      let yl = RealArray.length y in
      if n <> yl || n <> xl then raise IncompatibleArguments
    end;
    c_matvec rawptr x y

  let set_to_zero { payload = { data }; rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    Bigarray.Array2.fill data 0.0

  external c_copy : cptr -> t -> unit
    = "sunml_matrix_band_copy"

  let blit { rawptr = cptr1; valid = valid1 }
           ({ valid = valid2 } as m2) =
    if check_valid && not (valid1 && valid2) then raise Invalidated;
    c_copy cptr1 m2

  let space { payload = { dims = { n; mu; smu; ml } }; valid } =
    if check_valid && not valid then raise Invalidated;
    (n * (smu + ml + 1), 7 + n)

  let clone { payload = { data = dataa; dims }; valid } =
    if check_valid && not valid then raise Invalidated;
    create dims

  let ops = {
    m_clone      = clone;

    m_zero       = set_to_zero;

    m_copy       = blit;

    m_scale_add  = scale_add;

    m_scale_addi = scale_addi;

    m_matvec     = (fun { rawptr; valid } x y ->
                      if check_valid && not valid then raise Invalidated;
                      c_matvec rawptr x y);

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
      if m <= 0 then invalid_arg "m";
      if n <= 0 then invalid_arg "n";
      if nnz < 0 then invalid_arg "nnz"
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

  let unwrap ({ payload; rawptr } as m) =
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

  let set_col { payload = { idxptrs }; rawptr; valid } j idx =
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_set_idx rawptr j idx
    else idxptrs.{j} <- Index.of_int idx

  let get_col { payload = { idxptrs }; rawptr; valid } j =
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

  let set { payload = { idxvals; data }; rawptr; valid } idx i v =
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then (c_set_val rawptr idx i; c_set_data rawptr idx v)
    else (idxvals.{idx} <- Index.of_int i; data.{idx} <- v)

  let set_data { payload = { data }; rawptr; valid } idx v =
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_set_data rawptr idx v
    else data.{idx} <- v

  let set_rowval { payload = { idxvals }; rawptr; valid } idx i =
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_set_val rawptr idx i
    else idxvals.{idx} <- Index.of_int i

  let set_colval = set_rowval

  let get { payload = { idxvals; data }; rawptr; valid } idx =
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_get_val rawptr idx, c_get_data rawptr idx
    else Index.to_int idxvals.{idx}, data.{idx}

  let get_rowval { payload = { idxvals }; rawptr; valid } idx =
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_get_val rawptr idx
    else Index.to_int idxvals.{idx}

  let get_colval = get_rowval

  let get_data { payload = { data }; rawptr; valid } idx =
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_get_data rawptr idx
    else data.{idx}

  let pp (type s) fmt (mat : s t) =
    if check_valid && not mat.valid then raise Invalidated;
    let m, n = size mat in
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

  let ppi ?(start="[") ?(stop="]") ?(sep=";") ?(indent=4) ?(itemsep=" ")
          ?(rowcol=fun f->Format.fprintf f "%2d: ")
          ?(item=fun f->Format.fprintf f "%2d=% -15e") ()
          fmt mat =
    if check_valid && not mat.valid then raise Invalidated;
    let m, n = size mat in
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
        Format.pp_open_hovbox fmt 0;
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

  let scale_add c ({ rawptr = ptr1; valid = valid1 } as m1)
                   { rawptr = ptr2; valid = valid2 } =
    if check_valid && not (valid1 && valid2) then raise Invalidated;
    if Sundials_configuration.safe && c_size ptr1 <> c_size ptr2
    then raise IncompatibleArguments;
    c_scale_add c m1 ptr2

  external c_scale_addi : float -> 's t -> unit
    = "sunml_matrix_sparse_scale_addi"

  let scale_addi c ({ valid } as a) =
    if check_valid && not valid then raise Invalidated;
    c_scale_addi c a

  external c_matvec
    : cptr -> RealArray.t -> RealArray.t -> unit
    = "sunml_matrix_sparse_matvec"

  let matvec { rawptr; valid } x y =
    if check_valid && not valid then raise Invalidated;
    if Sundials_configuration.safe then begin
      let xl = RealArray.length x in
      let yl = RealArray.length y in
      let m, n = c_size rawptr in
      if n <> xl || m <> yl then raise IncompatibleArguments
    end;
    c_matvec rawptr x y

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

  let blit { rawptr = ptr1; valid = valid1 }
           ({ rawptr = ptr2; valid = valid2 } as m2) =
    if check_valid && not (valid1 && valid2) then raise Invalidated;
    if Sundials_configuration.safe && c_size ptr1 <> c_size ptr2
    then raise IncompatibleArguments;
    c_copy ptr1 m2

  external c_space : cptr -> int * int
      = "sunml_matrix_sparse_space"

  let space { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_space rawptr

  let clone { rawptr; valid; payload = {
      idxvals = valsa; idxptrs = ptrsa; data = dataa; sformat = fmta; } } =
    if check_valid && not valid then raise Invalidated;
    let m, n = c_size rawptr in
    let nnz, _ = c_dims rawptr in
    c_create m n nnz fmta

  let ops = {
    m_clone      = clone;

    m_zero       = set_to_zero;

    m_copy       = blit;

    m_scale_add  = scale_add;

    m_scale_addi = scale_addi;

    m_matvec     = (fun { rawptr; valid } x y ->
                      if check_valid && not valid then raise Invalidated;
                      c_matvec rawptr x y);

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
    let open RealArray2 in
    let m, n = size a in
    create m n

  let scale_addi c a =
    let ad = RealArray2.unwrap a in
    let an, am = Bigarray.Array2.(dim1 ad, dim2 ad) in
    for j = 0 to am - 1 do
      for i = 0 to an - 1 do
        ad.{i, j} <- c *. ad.{i, j} +. (if i = j then 1. else 0.)
      done
    done

  let scale_add c a b =
    let ad, bd = RealArray2.(unwrap a, unwrap b) in
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

    m_copy       = blit;

    m_scale_add  = scale_add;

    m_scale_addi = scale_addi;

    m_matvec     = matvec;

    m_space      = space;
  }
end (* }}} *)

module ArrayBand = struct (* {{{ *)

  type smu = int
  type mu = int
  type ml = int

  type t = RealArray2.t * (smu * mu * ml)

  let make ((smu, mu, ml) as dims) n v =
    (RealArray2.make (smu + ml + 1) n v, dims)

  let create ((smu, mu, ml) as dims) n =
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

  let unwrap (a, dims) = RealArray2.unwrap a

  external copy'
    : RealArray2.t -> RealArray2.t
      -> int * int * int * int -> unit
    = "sunml_arraybandmatrix_copy"

  let set_to_zero (x, dims) =
    Bigarray.Array2.fill (RealArray2.unwrap x) 0.0

  let blit (a, (a_smu, a_mu, a_ml)) (b, (b_smu, b_mu, b_ml))
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
    (RealArray2.create m n, dims)

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
    let ad, bd = RealArray2.(unwrap a, unwrap b) in
    let (an, am), (bn, bm) = RealArray2.(size a, size b) in
    if Sundials_configuration.safe && (am <> bm || mu_b > mu_a || ml_b > ml_a)
    then raise IncompatibleArguments;
    for j = 0 to bm - 1 do
      for i = - mu_b to ml_b do
        ad.{j, i + smu_a} <- c *. ad.{j, i + smu_a} +. bd.{j, i + smu_b}
      done;
    done

  let space (a, _) =
    let m, n = RealArray2.size a in
    (m * n, 5)
    (* 3 integer variables + 2 integer dimensions *)

  let ops = {
    m_clone      = clone;

    m_zero       = set_to_zero;

    m_copy       = blit;

    m_scale_add  = scale_add;

    m_scale_addi = scale_addi;

    m_matvec     = matvec;

    m_space      = space;
  }

end (* }}} *)

type standard
type custom

(* Must correspond with sundials_matrix_ml.h:mat_matrix_id_tag *)
type id =
  | Dense
  | Band
  | Sparse
  | Custom

type cmat

(* Must correspond with sundials_matrix_ml.h:mat_matrix_index *)
type ('k, 'm, 'nd, 'nk) t = {
  payload : 'm;
  rawptr  : cmat;
  id      : id;
  mat_ops : ('m, 'nd) matrix_ops;
}

type 'nk dense =
  (standard, Dense.t, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

type 'nk band =
  (standard, Band.t, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

type ('s, 'nk) sparse =
  (standard, 's Sparse.t, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

type 'nk arraydense = (custom, ArrayDense.t, RealArray.t, 'nk) t

type 'nk arrayband = (custom, ArrayBand.t, RealArray.t, 'nk) t

external c_wrap : id -> 'content_cptr -> 'm -> cmat
  = "sunml_matrix_wrap"

let wrap_dense (data : Dense.t) = {
    payload = data;
    rawptr  = c_wrap Dense data.rawptr data;
    id      = Dense;
    mat_ops = Dense.ops;
  }

let dense ?m ?(i=0.0) n =
  let m = match m with Some m -> m | None -> n in
  wrap_dense (Dense.make m n i)

let wrap_band (data : Band.t) = {
    payload = data;
    rawptr  = c_wrap Band data.rawptr data;
    id      = Band;
    mat_ops = Band.ops;
  }

let band ?(mu=2) ?smu ?ml ?(i=0.0) n =
  let ml  = match ml with Some ml -> ml | None -> mu in
  let smu = match smu with Some smu -> smu | None -> mu+ml in
  wrap_band (Band.(make { n; mu; smu; ml } i))

let wrap_sparse (data : 'f Sparse.t) =
  (match Config.sundials_version with
   | (2,v,_) when v < 6 -> raise Config.NotImplementedBySundialsVersion
   | _ -> ());
  {
    payload = data;
    rawptr  = c_wrap Sparse data.rawptr data;
    id      = Sparse;
    mat_ops = Sparse.ops;
  }

let sparse_csc ?m ?nnz n =
  let m = match m with Some m -> m | None -> n in
  let nnz = match nnz with Some nnz -> nnz | None -> n / 10 in
  wrap_sparse Sparse.(make CSC m n nnz)

let sparse_csr ?m ?nnz n =
  let m = match m with Some m -> m | None -> n in
  let nnz = match nnz with Some nnz -> nnz | None -> n / 10 in
  wrap_sparse Sparse.(make CSR m n nnz)

let wrap_custom ops data = {
    payload = data;
    rawptr  = c_wrap Custom ops data;
    id      = Custom;
    mat_ops = ops;
  }

let wrap_arraydense = wrap_custom ArrayDense.ops

let arraydense ?m ?(i=0.0) n =
  let m = match m with Some m -> m | None -> n in
  wrap_arraydense (ArrayDense.make m n i)

let wrap_arrayband = wrap_custom ArrayBand.ops

let arrayband ?mu ?smu ?ml ?(i=0.0) n =
  let mu = match mu, smu with
           | None, None  -> 2
           | Some mu, _  -> mu
           | _, Some smu -> smu
  in
  let smu = match smu with Some smu -> smu | None -> mu in
  let ml  = match ml  with Some ml  -> ml  | None -> mu in
  wrap_custom ArrayBand.ops (ArrayBand.make (smu, mu, ml) n i)

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

external c_matvec
  : ('k, 'm, 'nd, 'nk) t -> ('d, 'k) Nvector.t -> ('d, 'k) Nvector.t -> unit
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

let blit ({ payload = a; mat_ops = { m_copy } } as m1)
         ({ payload = b } as m2) =
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

(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "sunml_mat_init_module"

let _ =
  c_init_module
    (* Exceptions must be listed in the same order as
       mat_exn_index.  *)
    [|Invalidated;
      IncompatibleArguments;
      ZeroDiagonalElement 0;
    |]

