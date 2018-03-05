
exception Invalidated

let check_valid =
  match Sundials.sundials_version with
  | 2,_,_ -> Sundials_config.safe
  | _ -> false

(* Must correspond with matrix_ml.h:matrix_index *)
type ('data, 'cptr) matrix_content = {
  payload : 'data;
  rawptr  : 'cptr;
  mutable valid : bool;
}

type 'k serial_nvector
  = (Nvector_serial.data, [>Nvector_serial.kind] as 'k) Nvector.t

type ('m, 'd, 'k) matrix_ops = {
  m_clone      : 'm -> 'm;

  m_zero       : 'm -> unit;

  m_copy       : 'm -> 'm -> unit;

  m_scale_add  : float -> 'm -> 'm -> unit;

  m_scale_addi : float -> 'm -> unit;

  m_matvec     : 'm -> ('d, 'k) Nvector.t -> ('d, 'k) Nvector.t -> unit;

  m_space      : 'm -> int * int;
}

module Dense = struct

  type data = Sundials.RealArray2.data
  type cptr

  type t = (data, cptr) matrix_content

  external c_create : int -> int -> t
    = "c_matrix_dense_create"

  let create i j =
    if Sundials_config.safe && (i <= 0 || j <= 0)
    then failwith "Both M and N must be positive";
    c_create i j

  let unwrap { payload } = payload

  let make m n v =
    let { payload } as r = create m n in
    Bigarray.Array2.fill payload 0.0;
    r

  let invalidate v =
    if check_valid then v.valid <- false

  external c_size : cptr -> int * int
    = "c_matrix_dense_size"

  let size { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_size rawptr

  external c_print : Sundials.Logfile.t -> cptr -> unit
    = "c_matrix_dense_print"

  let print logfile { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_print logfile rawptr

  let pp fmt { payload=d; valid } =
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

  let ppi ?(start="[") ?(stop="]") ?(rowsep=";") ?(indent=4) ?(sep=" ")
          ?(item=fun f->Format.fprintf f "(%2d,%2d)=% -15e")
          fmt { payload=d; valid } =
    if check_valid && not valid then raise Invalidated;
    let ni, nj = Bigarray.Array2.dim2 d - 1, Bigarray.Array2.dim1 d - 1 in
    Format.pp_print_string fmt start;
    Format.pp_open_vbox fmt 0;
    for i = 0 to ni do

      Format.pp_open_hovbox fmt indent;
      for j = 0 to nj do
        if j > 0 then (
          Format.pp_print_string fmt sep;
          Format.pp_print_cut fmt ();
        );
        item fmt i j d.{j, i}
      done;
      Format.pp_close_box fmt ();

      if i < ni then (
        Format.pp_print_string fmt rowsep;
        Format.pp_print_cut fmt ();
      );

    done;
    Format.pp_close_box fmt ();
    Format.pp_print_string fmt stop

  (*
  external c_get : cptr -> int -> int -> float
    = "c_densematrix_get"

  let get { rawptr; valid } i j =
    if check_valid && not valid then raise Invalidated;
    c_get rawptr i j
  *)

  let get { payload; valid } i j =
    if check_valid && not valid then raise Invalidated;
    payload.{j, i}

  (*
  external c_set : cptr -> int -> int -> float -> unit
    = "c_densematrix_set"

  let set { rawptr; valid } i j e =
    if check_valid && not valid then raise Invalidated;
    c_set rawptr i j e
  *)

  let set { payload; valid } i j v =
    if check_valid && not valid then raise Invalidated;
    payload.{j, i} <- v

  let update { payload; valid } i j f =
    if check_valid && not valid then raise Invalidated;
    payload.{j, i} <- f payload.{j, i}

  external c_scale_add : float -> cptr -> cptr -> unit
    = "c_matrix_dense_scale_add"

  let scale_add c { rawptr = ptr1; valid = valid1 }
                  { rawptr = ptr2; valid = valid2 } =
    if check_valid && not (valid1 && valid2) then raise Invalidated;
    c_scale_add c ptr1 ptr2

  external c_scale_addi : float -> cptr -> unit
    = "c_matrix_dense_scale_addi"

  let scale_addi c { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_scale_addi c rawptr

  external c_matvec : cptr -> 'k serial_nvector -> 'k serial_nvector -> unit
    = "c_matrix_dense_matvec"

  let matvec { rawptr; valid } ~x ~y =
    if check_valid && not valid then raise Invalidated;
    c_matvec rawptr x y

  external c_zero : cptr -> unit
    = "c_matrix_dense_zero"

  let set_to_zero { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_zero rawptr

  external c_copy : cptr -> cptr -> unit
    = "c_matrix_dense_copy"

  let blit { rawptr = ptr1; valid = valid1 }
           { rawptr = ptr2; valid = valid2 } =
    if check_valid && not (valid1 && valid2) then raise Invalidated;
    c_copy ptr1 ptr2

  external c_space : cptr -> int * int
      = "c_matrix_dense_space"

  let space { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_space rawptr

  let clone { valid; payload = dataa } =
    if check_valid && not valid then raise Invalidated;
    let m = Bigarray.Array2.dim2 dataa in
    let n = Bigarray.Array2.dim1 dataa in
    let { payload = datab } as r = create m n in
    Bigarray.Array2.blit dataa datab;
    r

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
end

module Band = struct

  type data = {
    ismu : int;
    data : Sundials.RealArray2.data;
  }
  type cptr

  type t = (data, cptr) matrix_content

  type dimensions = {
      n   : int;
      mu  : int;
      smu : int;
      ml  : int;
    }

  external c_create : dimensions -> t
    = "c_matrix_band_create"

  let create ({n; mu; smu; ml} as dimensions) =
    if Sundials_config.safe && (n <= 0 || smu < 0 || ml < 0)
    then failwith "Dimensions must be positive";
    c_create dimensions

  let unwrap { payload = { data } } = data

  let make dimensions x =
    let { payload = { data } } as r = create dimensions in
    Bigarray.Array2.fill data x;
    r

  let invalidate v =
    if check_valid then v.valid <- false

  external c_size : cptr -> int * int
    = "c_matrix_band_size"

  let size { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_size rawptr

  external c_dims : cptr -> dimensions
      = "c_matrix_band_dims"

  let dims { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_dims rawptr

  external c_print : Sundials.Logfile.t -> cptr -> unit
    = "c_matrix_band_print"

  let print logfile { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_print logfile rawptr

  let pp fmt { payload={data=d}; valid } =
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

  let ppi ?(start="[") ?(stop="]") ?(rowsep=";") ?(indent=4) ?(sep=" ")
          ?(item=fun f->Format.fprintf f "(%2d,%2d)=% -15e")
          fmt { payload={data=d}; valid } =
    if check_valid && not valid then raise Invalidated;
    let ni, nj = Bigarray.Array2.dim2 d - 1, Bigarray.Array2.dim1 d - 1 in
    Format.pp_print_string fmt start;
    Format.pp_open_vbox fmt 0;
    for i = 0 to ni do

      Format.pp_open_hovbox fmt indent;
      for j = 0 to nj do
        if j > 0 then (
          Format.pp_print_string fmt sep;
          Format.pp_print_cut fmt ();
        );
        item fmt i j d.{j, i}
      done;
      Format.pp_close_box fmt ();

      if i < ni then (
        Format.pp_print_string fmt rowsep;
        Format.pp_print_cut fmt ();
      );

    done;
    Format.pp_close_box fmt ();
    Format.pp_print_string fmt stop

  (*
  external c_get : cptr -> int -> int -> float
    = "c_bandmatrix_get"

  let get { rawptr; valid } i j =
    if check_valid && not valid then raise Invalidated;
    c_get rawptr i j
  *)

  let get { payload = {data; ismu}; valid } i j =
    if check_valid && not valid then raise Invalidated;
    data.{j, i - j + ismu}

  (*
  external c_set : cptr -> int -> int -> float -> unit
    = "c_bandmatrix_set"

  let set { rawptr; valid } i j e =
    if check_valid && not valid then raise Invalidated;
    c_set rawptr i j e
  *)

  let set { payload; valid } i j v =
    if check_valid && not valid then raise Invalidated;
    payload.{j, i} <- v

  let update { payload; valid } i j f =
    if check_valid && not valid then raise Invalidated;
    payload.{j, i} <- f payload.{j, i}

  let set { payload={data; ismu}; valid } i j v =
    if check_valid && not valid then raise Invalidated;
    data.{j, i - j + ismu} <- v

  let update { payload={data; ismu}; valid } i j f =
    if check_valid && not valid then raise Invalidated;
    let k = i - j + ismu in
    data.{j, k} <- f data.{j, k}

  external c_scale_add : float -> cptr -> cptr -> unit
    = "c_matrix_band_scale_add"

  let scale_add c { rawptr = ptr1; valid = valid1 }
                  { rawptr = ptr2; valid = valid2 } =
    if check_valid && not (valid1 && valid2) then raise Invalidated;
    c_scale_add c ptr1 ptr2

  external c_scale_addi : float -> cptr -> unit
    = "c_matrix_band_scale_addi"

  let scale_addi c { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_scale_addi c rawptr

  external c_matvec : cptr -> 'a serial_nvector -> 'a serial_nvector -> unit
    = "c_matrix_band_matvec"

  let matvec { rawptr; valid } ~x ~y =
    if check_valid && not valid then raise Invalidated;
    c_matvec rawptr x y

  external c_zero : cptr -> unit
    = "c_matrix_band_zero"

  let set_to_zero { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_zero rawptr

  external c_copy : cptr -> cptr -> unit
    = "c_matrix_band_copy"

  let blit { rawptr = ptr1; valid = valid1 }
           { rawptr = ptr2; valid = valid2 } =
    if check_valid && not (valid1 && valid2) then raise Invalidated;
    c_copy ptr1 ptr2

  external c_space : cptr -> int * int
      = "c_matrix_band_space"

  let space { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_space rawptr

  let clone { rawptr; valid; payload = { data = dataa } } =
    if check_valid && not valid then raise Invalidated;
    let d = c_dims rawptr in
    let { payload = { data = datab } } as r = create d in
    Bigarray.Array2.blit dataa datab;
    r

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
end

module Sparse = struct
  type csc
  type csr

  type index_array =
    (int, Bigarray.int_elt, Bigarray.c_layout) Bigarray.Array1.t

  type sformat =
    | CSC_MAT
    | CSR_MAT

  type data = {
    idxvals : index_array;
    idxptrs : index_array;
    data    : Sundials.RealArray.t;
    sformat : sformat;
  }
    
  type cptr

  type 's t = (data, cptr) matrix_content

  type t_csc = csc t
  type t_csr = csr t

  external c_create : int -> int -> int -> sformat -> 'f t
    = "c_matrix_sparse_create"

  let make_csc m n nnz =
    if Sundials_config.safe && (m <= 0 || n <= 0 || nnz <= 0)
    then failwith "M, N and NNZ must all be positive";
    c_create m n nnz CSC_MAT

  let make_csr m n nnz =
    if Sundials_config.safe && (m <= 0 || n <= 0 || nnz <= 0)
    then failwith "M, N and NNZ must all be positive";
    c_create m n nnz CSR_MAT

  external c_from_dense : sformat -> Dense.cptr -> float -> 'f t
    = "c_matrix_sparse_from_dense"

  let csc_from_dense Dense.({ rawptr; valid }) droptol =
    if check_valid && not valid then raise Invalidated;
    c_from_dense CSC_MAT rawptr droptol

  let csr_from_dense Dense.({ rawptr; valid }) droptol =
    if check_valid && not valid then raise Invalidated;
    c_from_dense CSR_MAT rawptr droptol

  external c_from_band : sformat -> Band.cptr -> float -> 'f t
    = "c_matrix_sparse_from_band"

  let csc_from_band Band.({ rawptr; valid }) droptol =
    if check_valid && not valid then raise Invalidated;
    c_from_band CSC_MAT rawptr droptol

  let csr_from_band Band.({ rawptr; valid }) droptol =
    if check_valid && not valid then raise Invalidated;
    c_from_band CSR_MAT rawptr droptol

  external c_realloc : cptr -> unit
    = "c_matrix_sparse_realloc"

  let realloc { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_realloc rawptr

  let unwrap { payload = { idxvals; idxptrs; data } } = idxvals, idxptrs, data

  let invalidate v =
    if check_valid then v.valid <- false

  external c_size : cptr -> int * int
    = "c_matrix_sparse_size"

  let size { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_size rawptr

  external c_dims : cptr -> int * int
    = "c_matrix_sparse_dims"

  let dims { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_dims rawptr

  let set_col { payload = { idxptrs }; valid } j idx =
    if check_valid && not valid then raise Invalidated;
    idxptrs.{j} <- idx

  let get_col { payload = { idxptrs }; valid } j =
    if check_valid && not valid then raise Invalidated;
    idxptrs.{j}

  let set_row = set_col
  let get_row = get_col

  let set { payload = { idxvals; data }; valid } idx i v =
    if check_valid && not valid then raise Invalidated;
    idxvals.{idx} <- i;
    data.{idx} <- v

  let set_data { payload = { data }; valid } idx v =
    if check_valid && not valid then raise Invalidated;
    data.{idx} <- v

  let set_rowval { payload = { idxvals }; valid } idx i =
    if check_valid && not valid then raise Invalidated;
    idxvals.{idx} <- i

  let set_colval = set_rowval

  let get { payload = { idxvals; data }; valid } idx =
    if check_valid && not valid then raise Invalidated;
    idxvals.{idx}, data.{idx}

  let get_rowval { payload = { idxvals }; valid } idx =
    if check_valid && not valid then raise Invalidated;
    idxvals.{idx}

  let get_colval = get_rowval

  let get_data { payload = { data }; valid } idx =
    if check_valid && not valid then raise Invalidated;
    data.{idx}

  external c_print : Sundials.Logfile.t -> cptr -> unit
    = "c_matrix_sparse_print"

  let print logfile { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_print logfile rawptr

  let pp fmt mat =
    if check_valid && not mat.valid then raise Invalidated;
    let m, n = size mat in
    let end_row = ref false in
    let idxrole = match mat.payload.sformat with
                  | CSC_MAT -> "col"
                  | CSR_MAT -> "row" in

    Format.pp_print_string fmt "[";
    Format.pp_open_vbox fmt 0;
    for j = 0 to n - 1 do
      let p, np = get_col mat j, get_col mat (j + 1) in
      if p < np then begin
        if !end_row then begin
          Format.pp_print_string fmt ";";
          Format.pp_print_cut fmt ();
          Format.pp_close_box fmt ()
        end;
        Format.pp_open_hovbox fmt 4;
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

  let idxrole = function CSC_MAT -> "col" | CSR_MAT -> "row"

  let ppi ?(start="[") ?(stop="]") ?(rowsep=";") ?(indent=4) ?(sep=" ")
          ?(rowcol=fun f->Format.fprintf f "%2d: ")
          ?(item=fun f->Format.fprintf f "%2d=% -15e")
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
          Format.pp_print_string fmt rowsep;
          Format.pp_print_cut fmt ();
          Format.pp_close_box fmt ()
        end;
        Format.pp_open_hovbox fmt 4;
        rowcol fmt j;
        Format.pp_print_space fmt ();

        for i = p to np - 1 do
          if i > p then (
            Format.pp_print_string fmt sep;
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

  external c_scale_add : float -> cptr -> cptr -> unit
    = "c_matrix_sparse_scale_add"

  let scale_add c { rawptr = ptr1; valid = valid1 }
                  { rawptr = ptr2; valid = valid2 } =
    if check_valid && not (valid1 && valid2) then raise Invalidated;
    c_scale_add c ptr1 ptr2

  external c_scale_addi : float -> cptr -> unit
    = "c_matrix_sparse_scale_addi"

  let scale_addi c { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_scale_addi c rawptr

  external c_matvec : cptr -> 'a serial_nvector -> 'a serial_nvector -> unit
    = "c_matrix_sparse_matvec"

  let matvec { rawptr; valid } ~x ~y =
    if check_valid && not valid then raise Invalidated;
    c_matvec rawptr x y

  external c_zero : cptr -> unit
    = "c_matrix_sparse_zero"

  let set_to_zero { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_zero rawptr

  external c_copy : cptr -> cptr -> unit
    = "c_matrix_sparse_copy"

  let blit { rawptr = ptr1; valid = valid1 }
           { rawptr = ptr2; valid = valid2 } =
    if check_valid && not (valid1 && valid2) then raise Invalidated;
    c_copy ptr1 ptr2

  external c_space : cptr -> int * int
      = "c_matrix_sparse_space"

  let space { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_space rawptr

  let clone { rawptr; valid; payload = {
      idxvals = valsa; idxptrs = ptrsa; data = dataa; sformat = fmta; } } =
    if check_valid && not valid then raise Invalidated;
    let m, n = c_size rawptr in
    let nnz, _ = c_dims rawptr in
    let { payload = { idxvals = valsb; idxptrs = ptrsb; data = datab } } as r =
      c_create m n nnz fmta in
    Bigarray.Array1.blit valsa valsb;
    Bigarray.Array1.blit ptrsa ptrsb;
    Sundials.RealArray.blit dataa datab;
    r

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
end

type standard
type custom

type csc = Sparse.csc
type csr = Sparse.csr

type id =
  | Dense
  | Band
  | Sparse
  | Custom

type cptr

type ('k, 'm, 'nd, 'nk) t = {
  payload : 'm;
  rawptr  : cptr;
  id      : id;
  mat_ops : ('m, 'nd, 'nk) matrix_ops;
}

type 'nk dense =
  (standard, Dense.t, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

type 'nk band =
  (standard, Band.t, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

type ('s, 'nk) sparse =
  (standard, 's Sparse.t, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

external c_wrap_dense : Dense.t -> cptr
  = "c_matrix_wrap_dense"

let wrap_dense data = {
    payload = data;
    rawptr  = c_wrap_dense data;
    id      = Dense;
    mat_ops = Dense.ops;
  }

let make_dense m n x = wrap_dense (Dense.make m n x)

external c_wrap_band : Band.t -> cptr
  = "c_matrix_wrap_band"

let wrap_band data = {
    payload = data;
    rawptr  = c_wrap_band data;
    id      = Band;
    mat_ops = Band.ops;
  }

let make_band dims x = wrap_band (Band.make dims x)

external c_wrap_sparse : 's Sparse.t -> cptr
  = "c_matrix_wrap_sparse"

let wrap_sparse data = {
    payload = data;
    rawptr  = c_wrap_sparse data;
    id      = Sparse;
    mat_ops = Sparse.ops;
  }

let make_sparse_csc m n nnz = wrap_sparse (Sparse.make_csc m n nnz)

let make_sparse_csr m n nnz = wrap_sparse (Sparse.make_csr m n nnz)

external c_wrap_custom : ('m, 'd, 'k) matrix_ops -> 'm -> cptr
  = "c_matrix_wrap_custom"

let wrap_custom ops data = {
    payload = data;
    rawptr  = c_wrap_custom ops data;
    id      = Custom;
    mat_ops = ops;
  }

let get_ops { mat_ops } = mat_ops

let get_id { id } = id

let unwrap { payload } = payload

external c_scale_add : float -> cptr -> cptr -> unit
  = "c_matrix_scale_add"

let scale_add c { rawptr = ptr1; payload = a; mat_ops = { m_scale_add } }
                { rawptr = ptr2; payload = b } =
  match Sundials.sundials_version with
  | 2,_,_ -> m_scale_add c a b
  | _ -> c_scale_add c ptr1 ptr2

external c_scale_addi : float -> cptr -> unit
  = "c_matrix_scale_addi"

let scale_addi c { rawptr = ptr; payload = a; mat_ops = { m_scale_addi } } =
  match Sundials.sundials_version with
  | 2,_,_ -> m_scale_addi c a
  | _ -> c_scale_addi c ptr

external c_matvec : cptr -> ('d, 'k) Nvector.t -> ('d, 'k) Nvector.t -> unit
  = "c_matrix_matvec"

let matvec { rawptr = ptr; payload = a; mat_ops = { m_matvec } } ~x ~y =
  match Sundials.sundials_version with
  | 2,_,_ -> m_matvec a x y
  | _ -> c_matvec ptr x y

external c_zero : cptr -> unit
  = "c_matrix_zero"

let set_to_zero { rawptr = ptr; payload = a; mat_ops = { m_zero } } =
  match Sundials.sundials_version with
  | 2,_,_ -> m_zero a
  | _ -> c_zero ptr

external c_copy : cptr -> cptr -> unit
  = "c_matrix_copy"

let blit { rawptr = ptr1; payload = a; mat_ops = { m_copy } }
         { rawptr = ptr2; payload = b } =
  match Sundials.sundials_version with
  | 2,_,_ -> m_copy a b
  | _ -> c_copy ptr1 ptr2

external c_space : cptr -> int * int
    = "c_matrix_space"

let space { rawptr = ptr; payload = a; mat_ops = { m_space } } =
  match Sundials.sundials_version with
  | 2,_,_ -> m_space a
  | _ -> c_space ptr

