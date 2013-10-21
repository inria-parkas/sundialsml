type 'a pp = ?prec:int -> Format.formatter -> 'a -> unit
and  'a show = ?prec:int -> 'a -> string
and  'a dump = ?prec:int -> Format.formatter -> 'a -> unit
and  'a display = ?prec:int -> 'a -> string
and  'a ppout = ?prec:int -> 'a -> unit
and  'a pperr = ?prec:int -> 'a -> unit

let read_write_invariance = ref false

let with_read_write_invariance f =
  if !read_write_invariance then f ()
  else
    try
      read_write_invariance := true;
      let ret = f () in
      read_write_invariance := false;
      ret
    with exn -> read_write_invariance := false; raise exn

let pp_enclose left right needed pp fmt x =
  if needed then Format.fprintf fmt left;
  pp fmt x;
  if needed then Format.fprintf fmt right
let pp_parens needed pp fmt x = pp_enclose "(" ")" needed pp fmt x
let pp_brackets needed pp fmt x = pp_enclose "[" "]" needed pp fmt x
let pp_braces needed pp fmt x = pp_enclose "{" "}" needed pp fmt x
let pp_array_brackets needed pp fmt x = pp_enclose "[|" "|]" needed pp fmt x
let pp_double_quotes needed pp fmt x = pp_enclose "\"" "\"" needed pp fmt x
let pp_quotes needed pp fmt x = pp_enclose "'" "'" needed pp fmt x

let show_enclose left right needed str =
  if needed then left ^ str ^ right
  else str
let show_parens needed str = show_enclose "(" ")" needed str
let show_brackets needed str = show_enclose "[" "]" needed str
let show_braces needed str = show_enclose "{" "}" needed str
let show_array_brackets needed str = show_enclose "[|" "|]" needed str
let show_double_quotes needed str = show_enclose "\"" "\"" needed str
let show_quotes needed str = show_enclose "'" "'" needed str

module Prec = struct
  let app = 9
  let plus = 6
  let minus = plus
  let times = 7
  let div = 7
end

(* pp -> others *)
let dump_of_pp (pp : 'a pp) ?(prec=0) fmt x =
  with_read_write_invariance (fun () -> pp ~prec fmt x)
let show_of_pp (pp : 'a pp) ?(prec=0) x =
  let buf = Buffer.create 100 in
  let bfmt = Format.formatter_of_buffer buf in
  pp ~prec bfmt x;
  Format.pp_print_flush bfmt ();
  Buffer.contents buf
let display_of_pp (pp : 'a pp) ?(prec=0) x =
  with_read_write_invariance (fun () -> show_of_pp pp ~prec x)
let ppout_of_pp (pp : 'a pp) ?(prec=0) x =
  print_string (show_of_pp pp ~prec x)
let pperr_of_pp (pp : 'a pp) ?(prec=0) x =
  prerr_string (show_of_pp pp ~prec x)

(* show -> others *)
let pp_of_show (show : 'a show) ?(prec=0) fmt x =
  Format.pp_print_string fmt (show ~prec x)
let dump_of_show (show : 'a show) ?(prec=0) fmt x =
  with_read_write_invariance (fun () -> pp_of_show show ~prec fmt x)
let display_of_show (show : 'a show) ?(prec=0) x =
  with_read_write_invariance (fun () -> show ~prec x)
let ppout_of_show (show : 'a show) ?(prec=0) x =
  pp_of_show show ~prec Format.std_formatter x
let pperr_of_show (show : 'a show) ?(prec=0) x =
  pp_of_show show ~prec Format.err_formatter x

let printers_of_show show =
  pp_of_show show, dump_of_show show, show, display_of_show show,
  ppout_of_show show, pperr_of_show show

let printers_of_pp pp =
  pp, dump_of_pp pp, show_of_pp pp, display_of_pp pp,
  ppout_of_pp pp, pperr_of_pp pp

let pp_int, dump_int, show_int, display_int, ppout_int, pperr_int =
  printers_of_show (fun ?(prec=0) x ->
    show_parens (x < 0 && prec > 0) (string_of_int x))

let pp_unit, dump_unit, show_unit, display_unit, ppout_unit, pperr_unit =
  printers_of_show (fun ?(prec=0) () -> "()")

let pp_bool, dump_bool, show_bool, display_bool, ppout_bool, pperr_bool =
  printers_of_show (fun ?(prec=0) b -> string_of_bool b)

let pp_string, dump_string, show_string, display_string,
    ppout_string, pperr_string =
  printers_of_show (fun ?(prec=0) s ->
      show_double_quotes true (String.escaped s))

let pp_string_noquote, dump_string_noquote, show_string_noquote,
  display_string_noquote, ppout_string_noquote, pperr_string_noquote =
  printers_of_show (fun ?(prec=0) s -> s)

let pp_char, dump_char, show_char, display_char, ppout_char, pperr_char =
  printers_of_show (fun ?(prec=0) c -> show_quotes true (Char.escaped c))

let pp_char_noquote, dump_char_noquote, show_char_noquote,
  display_char_noquote, ppout_char_noquote, pperr_char_noquote =
  printers_of_show (fun ?(prec=0) s -> String.make 1 s)

let pp_float, dump_float, show_float, display_float, ppout_float, pperr_float =
  printers_of_show (fun ?(prec=0) f ->
    match classify_float f with
    | FP_nan -> "nan"
    | FP_infinite when !read_write_invariance ->
      if f > 0. then "infinity"
      else "neg_infinity"
    | FP_infinite | FP_normal | FP_subnormal | FP_zero ->
      show_parens (f < 0. && prec > 0)
      (if !read_write_invariance then string_of_float f
       else Printf.sprintf "%g" f)
  )

let pp_const ?(paren_prec=Prec.app + 1) s ?(prec=0) fmt x =
  pp_parens (paren_prec <= prec)
    (fun _ _ -> Format.pp_print_string fmt s)
    fmt x
let dump_const ?(paren_prec=0) s = dump_of_pp (pp_const ~paren_prec s)
let show_const ?(paren_prec=0) s = show_of_pp (pp_const ~paren_prec s)
let display_const ?(paren_prec=0) s = display_of_pp (pp_const ~paren_prec s)
let ppout_const ?(paren_prec=0) s = ppout_of_pp (pp_const ~paren_prec s)
let pperr_const ?(paren_prec=0) s = pperr_of_pp (pp_const ~paren_prec s)

(* Darn value restriction.  *)
let show_fun ?(prec=0) _ = "<fun>"
let pp_fun ?(prec=0) fmt f = pp_of_show show_fun ~prec fmt f
let dump_fun ?(prec=0) fmt f = dump_of_show show_fun ~prec fmt f
let display_fun ?(prec=0) f = display_of_show show_fun ~prec f
let ppout_fun ?(prec=0) f = ppout_of_show show_fun ~prec f
let pperr_fun ?(prec=0) f = pperr_of_show show_fun ~prec f


let pp_block opening delim closing fmt pp_elems =
  Format.fprintf fmt opening;
  let open Fstream in
  (match Fstream.decons pp_elems with
   | None -> ()
   | Some (pp_elem, pp_elems) ->
     pp_elem fmt;
     Fstream.iter (fun pp -> Format.fprintf fmt delim; pp fmt) pp_elems);
  Format.fprintf fmt closing

let pp_seq opening delim closing fmt pp_elems =
  pp_block ("@[" ^^ opening) (delim ^^ "@ ") (closing ^^ "@]") fmt pp_elems

let pp_fixarg (pp : 'a pp) x fmt = pp fmt x

type 'a pp_needed = Required of 'a | Optional of 'a

let pp_record pp_fields ?(prec=0) fmt x =
  let pp_field (name, (f : 'a pp)) fmt =
    Format.fprintf fmt "%s =@[<h2>@ " name;
    f ~prec fmt x;
    Format.fprintf fmt "@]"
  in
  let process_field = function
    | Required f -> Some (pp_field f)
    | Optional f when !read_write_invariance -> Some (pp_field f)
    | Optional _ -> None
  in
  pp_block "@[{@[<hv>" ";@ " "@]}@]" fmt
    (Fstream.filter_map process_field (Fstream.of_array pp_fields))

let const k _ = k

let pp_list_ix (pp : int -> 'a pp) ?(prec=0) fmt xs =
  pp_seq "[@[<hov>" ";" "@]]" fmt
    (Fstream.mapi (fun i x -> pp_fixarg (pp i) x)
       (Fstream.of_list xs))
let dump_list_ix pp_elem = dump_of_pp (pp_list_ix pp_elem)
let show_list_ix pp_elem = show_of_pp (pp_list_ix pp_elem)
let display_list_ix pp_elem = display_of_pp (pp_list_ix pp_elem)
let ppout_list_ix pp_elem = ppout_of_pp (pp_list_ix pp_elem)
let pperr_list_ix pp_elem = pperr_of_pp (pp_list_ix pp_elem)

let pp_list pp_elem ?(prec=0) fmt xs =
  pp_list_ix (const pp_elem) ~prec fmt xs
let dump_list pp_elem = dump_of_pp (pp_list pp_elem)
let show_list pp_elem = show_of_pp (pp_list pp_elem)
let display_list pp_elem = display_of_pp (pp_list pp_elem)
let ppout_list pp_elem = ppout_of_pp (pp_list pp_elem)
let pperr_list pp_elem = pperr_of_pp (pp_list pp_elem)

let pp_hvlist_ix (pp : int -> 'a pp) ?(prec=0) fmt xs =
  pp_seq "[@[<hv>" ";" "@]]" fmt
    (Fstream.mapi (fun i x -> pp_fixarg (pp i) x)
       (Fstream.of_list xs))
let dump_hvlist_ix pp_elem = dump_of_pp (pp_hvlist_ix pp_elem)
let show_hvlist_ix pp_elem = show_of_pp (pp_hvlist_ix pp_elem)
let display_hvlist_ix pp_elem = display_of_pp (pp_hvlist_ix pp_elem)
let ppout_hvlist_ix pp_elem = ppout_of_pp (pp_hvlist_ix pp_elem)
let pperr_hvlist_ix pp_elem = pperr_of_pp (pp_hvlist_ix pp_elem)

let pp_hvlist pp_elem ?(prec=0) fmt xs =
  pp_hvlist_ix (const pp_elem) ~prec fmt xs
let dump_hvlist pp_elem = dump_of_pp (pp_hvlist pp_elem)
let show_hvlist pp_elem = show_of_pp (pp_hvlist pp_elem)
let display_hvlist pp_elem = display_of_pp (pp_hvlist pp_elem)
let ppout_hvlist pp_elem = ppout_of_pp (pp_hvlist pp_elem)
let pperr_hvlist pp_elem = pperr_of_pp (pp_hvlist pp_elem)

let pp_vlist_ix (pp : int -> 'a pp) ?(prec=0) fmt xs =
  pp_seq "[@[<v>" ";" "@]]" fmt
    (Fstream.mapi (fun i x -> pp_fixarg (pp i) x)
       (Fstream.of_list xs))
let dump_vlist_ix pp_elem = dump_of_pp (pp_vlist_ix pp_elem)
let show_vlist_ix pp_elem = show_of_pp (pp_vlist_ix pp_elem)
let display_vlist_ix pp_elem = display_of_pp (pp_vlist_ix pp_elem)
let ppout_vlist_ix pp_elem = ppout_of_pp (pp_vlist_ix pp_elem)
let pperr_vlist_ix pp_elem = pperr_of_pp (pp_vlist_ix pp_elem)

let pp_vlist pp_elem ?(prec=0) fmt xs =
  pp_hvlist_ix (const pp_elem) ~prec fmt xs
let dump_vlist pp_elem = dump_of_pp (pp_vlist pp_elem)
let show_vlist pp_elem = show_of_pp (pp_vlist pp_elem)
let display_vlist pp_elem = display_of_pp (pp_vlist pp_elem)
let ppout_vlist pp_elem = ppout_of_pp (pp_vlist pp_elem)
let pperr_vlist pp_elem = pperr_of_pp (pp_vlist pp_elem)

let pp_array_like length get opening closing (pp : 'a pp) ?(prec=0) fmt xs
  =
  pp_seq opening ";" closing fmt
    (Fstream.map
       (fun i -> pp_fixarg pp (get xs i))
       (Fstream.enum 0 (length xs - 1)))

let pp_array pp_elem ?(prec=0) fmt xs =
  pp_array_like Array.length Array.get "[|@[<hov>" "@]|]" pp_elem fmt xs
let dump_array pp_elem = dump_of_pp (pp_array pp_elem)
let show_array pp_elem = show_of_pp (pp_array pp_elem)
let display_array pp_elem = display_of_pp (pp_array pp_elem)
let ppout_array pp_elem = ppout_of_pp (pp_array pp_elem)
let pperr_array pp_elem = pperr_of_pp (pp_array pp_elem)

let format_of_bigarray_kind kind =
    if Obj.magic kind = Bigarray.float32 then
      ("Bigarray.float32" : ('a,'b,'c) format)
    else if Obj.magic kind = Bigarray.float64 then
      ("Bigarray.float64" : ('a,'b,'c) format)
    else if Obj.magic kind = Bigarray.complex32 then
      ("Bigarray.complex32" : ('a,'b,'c) format)
    else if Obj.magic kind = Bigarray.complex64 then
      ("Bigarray.complex64" : ('a,'b,'c) format)
    else if Obj.magic kind = Bigarray.int8_signed then
      ("Bigarray.int8_signed" : ('a,'b,'c) format)
    else if Obj.magic kind = Bigarray.int8_unsigned then
      ("Bigarray.int8_unsigned" : ('a,'b,'c) format)
    else if Obj.magic kind = Bigarray.int16_signed then
      ("Bigarray.int16_signed" : ('a,'b,'c) format)
    else if Obj.magic kind = Bigarray.int16_unsigned then
      ("Bigarray.int16_unsigned" : ('a,'b,'c) format)
    else if Obj.magic kind = Bigarray.int then
      ("Bigarray.int" : ('a,'b,'c) format)
    else if Obj.magic kind = Bigarray.int32 then
      ("Bigarray.int32" : ('a,'b,'c) format)
    else if Obj.magic kind = Bigarray.int64 then
      ("Bigarray.int64" : ('a,'b,'c) format)
    else if Obj.magic kind = Bigarray.nativeint then
      ("Bigarray.nativeint" : ('a,'b,'c) format)
    else if Obj.magic kind = Bigarray.char then
      ("Bigarray.char" : ('a,'b,'c) format)
    else invalid_arg "Pprint.pp_bigarray1: Unrecognized bigarray kind"

let format_of_bigarray_layout layout =
    if Obj.magic layout = Bigarray.c_layout then
      ("Bigarray.c_layout" : ('a,'b,'c) format)
    else if Obj.magic layout = Bigarray.fortran_layout then
      ("Bigarray.fortran_layout" : ('a,'b,'c) format)
    else invalid_arg "Pprint.pp_bigarray1: Unrecognized bigarray layout"

let pp__fail : 'a. string -> 'a pp = fun fcn ?prec _ _ ->
  failwith (fcn ^ " shouldn't actually be called")

let pp_bigarray_kind _ _ ?prec fmt kind =
  Format.fprintf fmt (format_of_bigarray_kind kind)

let pp_bigarray_layout _ ?prec fmt layout =
  Format.fprintf fmt (format_of_bigarray_layout layout)

let pp_bigarray_c_layout = pp_bigarray_layout (pp__fail "")
let pp_bigarray_fortran_layout = pp_bigarray_layout (pp__fail "")

let pp_bigarray_float32_elt = pp__fail "pp_bigarray_float32_elt"
let pp_bigarray_float64_elt = pp__fail "pp_bigarray_float64_elt"
let pp_bigarray_complex32_elt = pp__fail "pp_bigarray_complex32_elt"
let pp_bigarray_complex64_elt = pp__fail "pp_bigarray_complex64_elt"
let pp_bigarray_int8_signed_elt = pp__fail "pp_bigarray_int8_signed_elt"
let pp_bigarray_int8_unsigned_elt = pp__fail "pp_bigarray_int8_unsigned_elt"
let pp_bigarray_int16_signed_elt = pp__fail "pp_bigarray_int16_signed_elt"
let pp_bigarray_int16_unsigned_elt = pp__fail "pp_bigarray_int16_unsigned_elt"
let pp_bigarray_int_elt = pp__fail "pp_bigarray_int_elt"
let pp_bigarray_int32_elt = pp__fail "pp_bigarray_int32_elt"
let pp_bigarray_int32_elt = pp__fail "pp_bigarray_int32_elt"
let pp_bigarray_int64_elt = pp__fail "pp_bigarray_int64_elt"
let pp_bigarray_nativeint_elt = pp__fail "pp_bigarray_nativeint_elt"
let pp_bigarray_char_elt = pp__fail "pp_bigarray_char_elt"

let pp_bigarray1 pp_elem _ _ ?(prec=0) fmt xs =
  if !read_write_invariance
  then
    let kind = format_of_bigarray_kind (Bigarray.Array1.kind xs)
    and layout = format_of_bigarray_layout (Bigarray.Array1.layout xs) in
    pp_parens (prec >= Prec.app)
      (pp_array_like Bigarray.Array1.dim Bigarray.Array1.get
         ("Bigarray.Array1.of_array "^^kind^^" "^^layout^^" [|")
         "|]" pp_elem ~prec:Prec.app) fmt xs
  else
    pp_array_like Bigarray.Array1.dim Bigarray.Array1.get
      "[|" "|]" pp_elem ~prec fmt xs

let dump_bigarray1 pp_elem pp_repr pp_layout =
  dump_of_pp (pp_bigarray1 pp_elem pp_repr pp_layout)
let show_bigarray1 pp_elem pp_repr pp_layout =
  show_of_pp (pp_bigarray1 pp_elem pp_repr pp_layout)
let display_bigarray1 pp_elem pp_repr pp_layout =
  display_of_pp (pp_bigarray1 pp_elem pp_repr pp_layout)
let ppout_bigarray1 pp_elem pp_repr pp_layout =
  ppout_of_pp (pp_bigarray1 pp_elem pp_repr pp_layout)
let pperr_bigarray1 pp_elem pp_repr pp_layout =
  pperr_of_pp (pp_bigarray1 pp_elem pp_repr pp_layout)

let pp_tuple pp_fields ?(prec=0) fmt x =
  pp_seq "(@[<hv>" "," "@])" fmt
    (Fstream.map (fun pp -> pp_fixarg pp x)
       (Fstream.of_list pp_fields))
let show_tuple pp_fields = show_of_pp (pp_tuple pp_fields)
let dump_tuple pp_fields = dump_of_pp (pp_tuple pp_fields)
let display_tuple pp_fields = display_of_pp (pp_tuple pp_fields)
let ppout_tuple pp_fields = ppout_of_pp (pp_tuple pp_fields)
let pperr_tuple pp_fields = pperr_of_pp (pp_tuple pp_fields)

let pp_pair (pp1 : 'a pp) (pp2 : 'b pp) =
  pp_tuple [(fun ?prec fmt (x,y) -> pp1 ~prec:0 fmt x);
            (fun ?prec fmt (x,y) -> pp2 ~prec:0 fmt y)]
let show_pair pp1 pp2 = show_of_pp (pp_pair pp1 pp2)
let dump_pair pp1 pp2 = dump_of_pp (pp_pair pp1 pp2)
let display_pair pp1 pp2 = display_of_pp (pp_pair pp1 pp2)
let ppout_pair pp1 pp2 = ppout_of_pp (pp_pair pp1 pp2)
let pperr_pair pp1 pp2 = pperr_of_pp (pp_pair pp1 pp2)

let pp_triple (pp1: 'a pp) (pp2 : 'b pp) (pp3 : 'c pp) =
  pp_tuple [(fun ?prec fmt (x,y,z) -> pp1 fmt x);
            (fun ?prec fmt (x,y,z) -> pp2 fmt y);
            (fun ?prec fmt (x,y,z) -> pp3 fmt z)]
let show_triple pp1 pp2 pp3 = show_of_pp (pp_triple pp1 pp2 pp3)
let dump_triple pp1 pp2 pp3 = dump_of_pp (pp_triple pp1 pp2 pp3)
let display_triple pp1 pp2 pp3 = display_of_pp (pp_triple pp1 pp2 pp3)
let ppout_triple pp1 pp2 pp3 = ppout_of_pp (pp_triple pp1 pp2 pp3)
let pperr_triple pp1 pp2 pp3 = pperr_of_pp (pp_triple pp1 pp2 pp3)

let pp_precompose f (pp : 'a pp) ?prec fmt x = pp ?prec fmt (f x)
let show_precompose f pp = show_of_pp (pp_precompose f pp)
let dump_precompose f pp = dump_of_pp (pp_precompose f pp)
let display_precompose f pp = display_of_pp (pp_precompose f pp)
let ppout_precompose f pp = ppout_of_pp (pp_precompose f pp)
let pperr_precompose f pp = pperr_of_pp (pp_precompose f pp)

let pp_fst pp = pp_precompose fst pp
let show_fst pp = show_of_pp (pp_fst pp)
let dump_fst pp = dump_of_pp (pp_fst pp)
let display_fst pp = display_of_pp (pp_fst pp)
let ppout_fst pp = ppout_of_pp (pp_fst pp)
let pperr_fst pp = pperr_of_pp (pp_fst pp)

let pp_snd pp = pp_precompose snd pp
let show_snd pp = show_of_pp (pp_snd pp)
let dump_snd pp = dump_of_pp (pp_snd pp)
let display_snd pp = display_of_pp (pp_snd pp)
let ppout_snd pp = ppout_of_pp (pp_snd pp)
let pperr_snd pp = pperr_of_pp (pp_snd pp)

let pp_ctor name (pp : 'a pp) ?(prec=0) fmt x =
  match name with
  | Optional name when not !read_write_invariance -> pp ~prec fmt x
  | Required name | Optional name ->
    pp_parens (prec >= Prec.app)
      (fun fmt x ->
         Format.fprintf fmt "%s " name;
         pp ~prec:Prec.app fmt x)
      fmt x
let show_ctor name pp = show_of_pp (pp_ctor name pp)
let dump_ctor name pp = dump_of_pp (pp_ctor name pp)
let display_ctor name pp = display_of_pp (pp_ctor name pp)
let ppout_ctor name pp = ppout_of_pp (pp_ctor name pp)
let pperr_ctor name pp = pperr_of_pp (pp_ctor name pp)

let pp_option pp_elem ?(prec=0) fmt = function
  | Some x -> pp_ctor (Required "Some") pp_elem ~prec fmt x
  | None -> pp_string_noquote fmt "None"
let show_option pp_elem = show_of_pp (pp_option pp_elem)
let dump_option pp_elem = dump_of_pp (pp_option pp_elem)
let display_option pp_elem = display_of_pp (pp_option pp_elem)
let ppout_option pp_elem = ppout_of_pp (pp_option pp_elem)
let pperr_option pp_elem = pperr_of_pp (pp_option pp_elem)

let pp_lazy_t pp_elem ?(prec=0) fmt x =
  if Lazy.lazy_is_val x
  then pp_ctor (Required "lazy") pp_elem ~prec fmt (Lazy.force x)
  else pp_string_noquote fmt "<lazy>"
let show_lazy_t pp_elem = show_of_pp (pp_lazy_t pp_elem)
let dump_lazy_t pp_elem = dump_of_pp (pp_lazy_t pp_elem)
let display_lazy_t pp_elem = display_of_pp (pp_lazy_t pp_elem)
let ppout_lazy_t pp_elem = ppout_of_pp (pp_lazy_t pp_elem)
let pperr_lazy_t pp_elem = pperr_of_pp (pp_lazy_t pp_elem)

let pp_lazy pp_elem = pp_lazy_t pp_elem
let show_lazy pp_elem = show_lazy_t pp_elem
let dump_lazy pp_elem = dump_lazy_t pp_elem
let display_lazy pp_elem = display_lazy_t pp_elem
let ppout_lazy pp_elem = ppout_lazy_t pp_elem
let pperr_lazy pp_elem = pperr_lazy_t pp_elem

let pp_exn, dump_exn, show_exn, display_exn, ppout_exn, pperr_exn =
  printers_of_show (fun ?(prec=0) e ->
      let str = Printexc.to_string e in
      let need_paren = prec >= Prec.app &&
                       (String.contains str ' ' || String.contains str '(')
      in show_parens need_paren str
    )
