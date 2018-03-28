
exception Invalidated
exception IncompatibleArguments

let check_valid =
  match Sundials.sundials_version with
  | 2,_,_ -> Sundials_config.safe
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

  type data = Sundials.real_array2
  type cptr

  type t = (data, cptr) matrix_content

  external c_create : int -> int -> t
    = "ml_matrix_dense_create"

  let create i j =
    if Sundials_config.safe then begin
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
          ?(item=fun f->Format.fprintf f "(%2d,%2d)=% -15e")
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
    = "ml_matrix_dense_scale_add"

  let scale_add c { rawptr = ptr1; valid = valid1 }
                  { rawptr = ptr2; valid = valid2 } =
    if check_valid && not (valid1 && valid2) then raise Invalidated;
    (* compatability checked on C side *)
    c_scale_add c ptr1 ptr2

  external c_scale_addi : float -> cptr -> unit
    = "ml_matrix_dense_scale_addi"

  let scale_addi c { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_scale_addi c rawptr

  external c_matvec : cptr -> 'k serial_nvector -> 'k serial_nvector -> unit
    = "ml_matrix_dense_matvec"

  let matvec { rawptr; payload; valid } ~x ~y =
    if check_valid && not valid then raise Invalidated;
    if Sundials_config.safe then begin
      let xl = Sundials.RealArray.length (Nvector.unwrap x) in
      let yl = Sundials.RealArray.length (Nvector.unwrap y) in
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
      = "ml_matrix_dense_space"

  let space { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_space rawptr

  let clone { payload = dataa; valid } =
    if check_valid && not valid then raise Invalidated;
    let m, n = Bigarray.Array2.(dim2 dataa, dim1 dataa) in
    let { payload = datab; rawptr = cptrb } as r = create m n in
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

  (* Must correspond with matrix_ml.h:mat_band_dimensions_index *)
  type dimensions = {
      n   : int;
      mu  : int;
      smu : int;
      ml  : int;
    }

  (* Must correspond with matrix_ml.h:mat_band_data_index *)
  type data = {
    data : Sundials.real_array2;
    dims : dimensions;
  }
  type cptr

  type t = (data, cptr) matrix_content

  external c_create : dimensions -> t
    = "ml_matrix_band_create"

  let create ({n; mu; smu; ml} as dimensions) =
    if Sundials_config.safe then begin
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

  let pp fmt { payload={data = d}; rawptr; valid } =
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
          ?(item=fun f->Format.fprintf f "(%2d,%2d)=% -15e")
          fmt { payload={data = d}; rawptr; valid } =
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
    = "ml_matrix_band_scale_add"
  let _ = Callback.register "ml_matrix_band_scale_add" c_scale_add

  let scale_add c ({ valid = valid1 } as m1)
                  { rawptr = ptr2; valid = valid2 } =
    if check_valid && not (valid1 && valid2) then raise Invalidated;
    c_scale_add c m1 ptr2

  external c_scale_addi : float -> cptr -> unit
    = "ml_matrix_band_scale_addi"

  let scale_addi c { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_scale_addi c rawptr

  external c_matvec : cptr -> 'a serial_nvector -> 'a serial_nvector -> unit
    = "ml_matrix_band_matvec"

  let matvec { rawptr; payload = { data; dims = { n } }; valid } ~x ~y =
    if check_valid && not valid then raise Invalidated;
    if Sundials_config.safe then begin
      let xl = Sundials.RealArray.length (Nvector.unwrap x) in
      let yl = Sundials.RealArray.length (Nvector.unwrap y) in
      if n <> yl || n <> xl then raise IncompatibleArguments
    end;
    c_matvec rawptr x y

  let set_to_zero { payload = { data }; rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    Bigarray.Array2.fill data 0.0

  external c_copy : cptr -> t -> unit
    = "ml_matrix_band_copy"
  let _ = Callback.register "ml_matrix_band_copy" c_copy

  let blit { rawptr = cptr1; valid = valid1 }
           ({ valid = valid2 } as m2) =
    if check_valid && not (valid1 && valid2) then raise Invalidated;
    c_copy cptr1 m2

  let space { payload = { dims = { n; mu; smu; ml } }; valid } =
    if check_valid && not valid then raise Invalidated;
    (n * (smu + ml + 1), 7 + n)

  let clone { payload = { data = dataa; dims }; valid } =
    if check_valid && not valid then raise Invalidated;
    let { payload = { data = datab } } as r = create dims in
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

  (* The payload field cannot be relied upon for Sundials < 3.0.0 *)
  let unsafe_content =
    match Sundials.sundials_version with
    | 2,_,_ -> true
    | _ -> false

  type index_array =
    (int, Bigarray.int_elt, Bigarray.c_layout) Bigarray.Array1.t

  (* See MAT_TO/FROM_SFORMAT(x) in matrix_ml.h *)
  type sformat =
    | CSC_MAT   (* Must correpond with Sundials CSC_MAT = 0 constant. *)
    | CSR_MAT   (* Must correpond with Sundials CSR_MAT = 1 constant. *)

  (* Must correspond with matrix_ml.h:mat_sparse_data_index *)
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
    = "ml_matrix_sparse_create"

  let make_csc m n nnz =
    if Sundials_config.safe then begin
      if m <= 0 then invalid_arg "m";
      if n <= 0 then invalid_arg "n";
      if nnz <= 0 then invalid_arg "nnz"
    end;
    c_create m n nnz CSC_MAT

  let make_csr m n nnz =
    if Sundials_config.safe then begin
      if m <= 0 then invalid_arg "m";
      if n <= 0 then invalid_arg "n";
      if nnz <= 0 then invalid_arg "nnz"
    end;
    if Sundials_config.safe then
      (match Sundials.sundials_version with
       | 2,v,_ when v < 7 -> raise Sundials.NotImplementedBySundialsVersion
       | _ -> ());
    c_create m n nnz CSR_MAT

  external c_from_dense : sformat -> Dense.cptr -> float -> 'f t
    = "ml_matrix_sparse_from_dense"

  let csc_from_dense droptol Dense.({ rawptr; valid }) =
    if check_valid && not valid then raise Invalidated;
    if Sundials_config.safe && droptol < 0.0 then invalid_arg "droptol";
    c_from_dense CSC_MAT rawptr droptol

  let csr_from_dense droptol Dense.({ rawptr; valid }) =
    if check_valid && not valid then raise Invalidated;
    if Sundials_config.safe && droptol < 0.0 then invalid_arg "droptol";
    if Sundials_config.safe then
      (match Sundials.sundials_version with
       | 2,v,_ when v < 7 -> raise Sundials.NotImplementedBySundialsVersion
       | _ -> ());
    c_from_dense CSR_MAT rawptr droptol

  external c_from_band : sformat -> Band.cptr -> float -> 'f t
    = "ml_matrix_sparse_from_band"

  let csc_from_band droptol Band.({ rawptr; valid }) =
    if check_valid && not valid then raise Invalidated;
    if Sundials_config.safe && droptol < 0.0 then invalid_arg "droptol";
    c_from_band CSC_MAT rawptr droptol

  let csr_from_band droptol Band.({ rawptr; valid }) =
    if check_valid && not valid then raise Invalidated;
    if Sundials_config.safe && droptol < 0.0 then invalid_arg "droptol";
    if Sundials_config.safe then
      (match Sundials.sundials_version with
       | 2,v,_ when v < 7 -> raise Sundials.NotImplementedBySundialsVersion
       | _ -> ());
    c_from_band CSR_MAT rawptr droptol

  external c_rewrap : 'f t -> data
    = "ml_matrix_sparse_rewrap"

  let unwrap ({ payload; rawptr } as m) =
    let { idxvals; idxptrs; data } =
      if unsafe_content then c_rewrap m else payload
    in
    idxvals, idxptrs, data

  external c_resize : 's t -> int -> bool -> unit
    = "ml_matrix_sparse_resize"

  let resize ?nnz a =
    let nnz = match nnz with Some x -> x | None -> 0 in
    c_resize a nnz true

  let invalidate v =
    if check_valid then v.valid <- false

  external c_size : cptr -> int * int
    = "ml_matrix_sparse_size"

  let size { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_size rawptr

  external c_dims : cptr -> int * int
    = "ml_matrix_sparse_dims"

  let dims { rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    c_dims rawptr

  external c_set_idx : cptr -> int -> int -> unit
    = "ml_matrix_sparse_set_idx"

  external c_get_idx : cptr -> int -> int
    = "ml_matrix_sparse_get_idx"

  let set_col { payload = { idxptrs }; rawptr; valid } j idx =
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_set_idx rawptr j idx
    else idxptrs.{j} <- idx

  let get_col { payload = { idxptrs }; rawptr; valid } j =
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_get_idx rawptr j
    else idxptrs.{j}

  let set_row = set_col
  let get_row = get_col

  external c_set_data : cptr -> int -> float -> unit
    = "ml_matrix_sparse_set_data"

  external c_get_data : cptr -> int -> float
    = "ml_matrix_sparse_get_data"

  external c_set_val : cptr -> int -> int -> unit
    = "ml_matrix_sparse_set_val"

  external c_get_val : cptr -> int -> int
    = "ml_matrix_sparse_get_val"

  let set { payload = { idxvals; data }; rawptr; valid } idx i v =
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then (c_set_val rawptr idx i; c_set_data rawptr idx v)
    else (idxvals.{idx} <- i; data.{idx} <- v)

  let set_data { payload = { data }; rawptr; valid } idx v =
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_set_data rawptr idx v
    else data.{idx} <- v

  let set_rowval { payload = { idxvals }; rawptr; valid } idx i =
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_set_val rawptr idx i
    else idxvals.{idx} <- i

  let set_colval = set_rowval

  let get { payload = { idxvals; data }; rawptr; valid } idx =
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_get_val rawptr idx, c_get_data rawptr idx
    else idxvals.{idx}, data.{idx}

  let get_rowval { payload = { idxvals }; rawptr; valid } idx =
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_get_val rawptr idx
    else idxvals.{idx}

  let get_colval = get_rowval

  let get_data { payload = { data }; rawptr; valid } idx =
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_get_data rawptr idx
    else data.{idx}

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

  let ppi ?(start="[") ?(stop="]") ?(sep=";") ?(indent=4) ?(itemsep=" ")
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
          Format.pp_print_string fmt sep;
          Format.pp_print_cut fmt ();
          Format.pp_close_box fmt ()
        end;
        Format.pp_open_hovbox fmt 4;
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
    = "ml_matrix_sparse_scale_add"
  let _ = Callback.register "ml_matrix_sparse_scale_add" c_scale_add

  let scale_add c ({ rawptr = ptr1; valid = valid1 } as m1)
                  { rawptr = ptr2; valid = valid2 } =
    if check_valid && not (valid1 && valid2) then raise Invalidated;
    if Sundials_config.safe && c_size ptr1 <> c_size ptr2
    then raise IncompatibleArguments;
    c_scale_add c m1 ptr2

  external c_scale_addi : float -> 's t -> unit
    = "ml_matrix_sparse_scale_addi"
  let _ = Callback.register "ml_matrix_sparse_scale_addi" c_scale_addi

  let scale_addi c ({ valid } as a) =
    if check_valid && not valid then raise Invalidated;
    c_scale_addi c a

  external c_matvec : cptr -> 'a serial_nvector -> 'a serial_nvector -> unit
    = "ml_matrix_sparse_matvec"

  let matvec { rawptr; valid } ~x ~y =
    if check_valid && not valid then raise Invalidated;
    if Sundials_config.safe then begin
      let xl = Sundials.RealArray.length (Nvector.unwrap x) in
      let yl = Sundials.RealArray.length (Nvector.unwrap y) in
      let m, n = c_size rawptr in
      if n <> xl || m <> yl then raise IncompatibleArguments
    end;
    c_matvec rawptr x y

  external c_set_to_zero : cptr -> unit
    = "ml_matrix_sparse_set_to_zero"

  let set_to_zero { payload = { idxvals; idxptrs; data }; rawptr; valid } =
    if check_valid && not valid then raise Invalidated;
    if unsafe_content then c_set_to_zero rawptr
    else (Bigarray.Array1.fill idxvals 0;
          Bigarray.Array1.fill idxptrs 0;
          Bigarray.Array1.fill data 0.0)

  external c_copy : cptr -> 's t -> unit
      = "ml_matrix_sparse_copy"
  let _ = Callback.register "ml_matrix_sparse_copy" c_copy

  let blit { rawptr = ptr1; valid = valid1 }
           ({ rawptr = ptr2; valid = valid2 } as m2) =
    if check_valid && not (valid1 && valid2) then raise Invalidated;
    if Sundials_config.safe && c_size ptr1 <> c_size ptr2
    then raise IncompatibleArguments;
    c_copy ptr1 m2

  external c_space : cptr -> int * int
      = "ml_matrix_sparse_space"

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

(* Must correspond with matrix_ml.h:mat_matrix_id_tag *)
type id =
  | Dense
  | Band
  | Sparse
  | Custom

type cptr

(* Must correspond with matrix_ml.h:mat_matrix_index *)
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

external c_wrap : id -> 'a -> cptr
  = "ml_matrix_wrap"

let wrap_dense data = {
    payload = data;
    rawptr  = c_wrap Dense data;
    id      = Dense;
    mat_ops = Dense.ops;
  }

let make_dense m n x = wrap_dense (Dense.make m n x)

let wrap_band data = {
    payload = data;
    rawptr  = c_wrap Band data;
    id      = Band;
    mat_ops = Band.ops;
  }

let make_band dims x = wrap_band (Band.make dims x)

let wrap_sparse data = {
    payload = data;
    rawptr  = c_wrap Sparse data;
    id      = Sparse;
    mat_ops = Sparse.ops;
  }

let make_sparse_csc m n nnz = wrap_sparse (Sparse.make_csc m n nnz)

let make_sparse_csr m n nnz = wrap_sparse (Sparse.make_csr m n nnz)

external c_wrap_custom : ('m, 'd, 'k) matrix_ops -> 'm -> cptr
  = "ml_matrix_wrap_custom"

let wrap_custom ops data = {
    payload = data;
    rawptr  = c_wrap_custom ops data;
    id      = Custom;
    mat_ops = ops;
  }

let get_ops { mat_ops } = mat_ops

let get_id { id } = id

let unwrap { payload } = payload

external c_scale_add
  : float -> ('k, 'm, 'nd, 'nk) t -> ('k, 'm, 'nd, 'nk) t -> unit
  = "ml_matrix_scale_add"

let scale_add c ({ payload = a; mat_ops = { m_scale_add } } as m1)
                ({ payload = b } as m2) =
  match Sundials.sundials_version with
  | 2,_,_ -> m_scale_add c a b
  | _ -> c_scale_add c m1 m2

external c_scale_addi : float -> ('k, 'm, 'nd, 'nk) t -> unit
  = "ml_matrix_scale_addi"

let scale_addi c ({ payload = a; mat_ops = { m_scale_addi } } as m) =
  match Sundials.sundials_version with
  | 2,_,_ -> m_scale_addi c a
  | _ -> c_scale_addi c m

external c_matvec
  : ('k, 'm, 'nd, 'nk) t -> ('d, 'k) Nvector.t -> ('d, 'k) Nvector.t -> unit
  = "ml_matrix_matvec"

let matvec ({ payload = a; mat_ops = { m_matvec } } as m) ~x ~y =
  match Sundials.sundials_version with
  | 2,_,_ -> m_matvec a x y
  | _ -> c_matvec m x y

external c_zero : ('k, 'm, 'nd, 'nk) t -> unit
  = "ml_matrix_zero"

let set_to_zero ({ payload = a; mat_ops = { m_zero } } as m) =
  match Sundials.sundials_version with
  | 2,_,_ -> m_zero a
  | _ -> c_zero m

external c_copy : ('k, 'm, 'nd, 'nk) t -> ('k, 'm, 'nd, 'nk) t -> unit
  = "ml_matrix_copy"

let blit ({ payload = a; mat_ops = { m_copy } } as m1)
         ({ payload = b } as m2) =
  match Sundials.sundials_version with
  | 2,_,_ -> m_copy a b
  | _ -> c_copy m1 m2

external c_space : ('k, 'm, 'nd, 'nk) t -> int * int
    = "ml_matrix_space"

let space ({ payload = a; mat_ops = { m_space } } as m) =
  match Sundials.sundials_version with
  | 2,_,_ -> m_space a
  | _ -> c_space m

(* Let C code know about some of the values in this module.  *)
external c_init_module : exn array -> unit =
  "ml_mat_init_module"

let _ =
  c_init_module
    (* Exceptions must be listed in the same order as
       mat_exn_index.  *)
    [|Invalidated;
      IncompatibleArguments;
    |]

