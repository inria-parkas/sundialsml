type 'a pp = Format.formatter -> 'a -> unit
and  'a show = 'a -> string
and  'a dump = 'a pp
and  'a display = 'a show
and  'a print = 'a -> unit
and  'a prerr = 'a -> unit

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

(* pp -> others *)
let dump_of_pp pp fmt x = with_read_write_invariance (fun () -> pp fmt x)
let show_of_pp pp x =
  let buf = Buffer.create 100 in
  let bfmt = Format.formatter_of_buffer buf in
  pp bfmt x;
  Format.pp_print_flush bfmt ();
  Buffer.contents buf
let display_of_pp pp x =
  with_read_write_invariance (fun () -> show_of_pp pp x)
let print_of_pp pp x =
  pp Format.std_formatter x;
  (* Flushing every time isn't very good performance-wise, but it seems to be
     necessary so that this function's outputs interleaves properly with direct
     outputs to the stdout channel.  *)
  Format.pp_print_flush Format.std_formatter ()
let prerr_of_pp pp x =
  pp Format.err_formatter x;
  (* Flushing every time isn't very good performance-wise, but it seems to be
     necessary so that this function's outputs interleaves properly with direct
     outputs to the stderr channel.  *)
  Format.pp_print_flush Format.std_formatter ()

(* show -> others *)
let pp_of_show show fmt x = Format.fprintf fmt "%s" (show x)
let dump_of_show show fmt x =
  with_read_write_invariance (fun () -> pp_of_show show fmt x)
let display_of_show show x =
  with_read_write_invariance (fun () -> show x)
let print_of_show show x =
  pp_of_show show Format.std_formatter x;
  (* Flushing every time isn't very good performance-wise, but it seems to be
     necessary so that this function's outputs interleaves properly with direct
     outputs to the stdout channel.  *)
  Format.pp_print_flush Format.std_formatter ()
let prerr_of_show show x =
  pp_of_show show Format.err_formatter x;
  (* Flushing every time isn't very good performance-wise, but it seems to be
     necessary so that this function's outputs interleaves properly with direct
     outputs to the stderr channel.  *)
  Format.pp_print_flush Format.err_formatter ()

let printers_of_show show =
  pp_of_show show, dump_of_show show, show, display_of_show show,
  print_of_show show, prerr_of_show show

let printers_of_pp pp =
  pp, dump_of_pp pp, show_of_pp pp, display_of_pp pp,
  print_of_pp pp, prerr_of_pp pp


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


let pp_int, dump_int, show_int, display_int, _, _ =
  printers_of_show string_of_int

let pp_bool, dump_bool, show_bool, display_bool, _, _ =
  printers_of_show string_of_bool

let pp_string, dump_string, show_string, display_string,
    print_string, prerr_string =
  printers_of_show (fun s -> show_double_quotes true (String.escaped s))

let pp_string_verbatim, dump_string_verbatim, show_string_verbatim,
  display_string_verbatim, print_string_verbatim, prerr_string_verbatim =
  printers_of_show (fun s -> s)

let pp_char, dump_char, show_char, display_char, print_char, prerr_char =
  printers_of_show (fun c -> show_quotes true (Char.escaped c))

let pp_char_verbatim, dump_char_verbatim, show_char_verbatim,
  display_char_verbatim, print_char_verbatim, prerr_char_verbatim =
  printers_of_show (fun s -> String.make 1 s)

let pp_float, dump_float, show_float, display_float, print_float, prerr_float =
  printers_of_show (fun f ->
    match classify_float f with
    | FP_nan -> "nan"
    | FP_infinite when !read_write_invariance ->
      if f > 0. then "infinity"
      else "neg_infinity"
    | FP_infinite | FP_normal | FP_subnormal | FP_zero ->
      if !read_write_invariance then string_of_float f
      else Printf.sprintf "%g" f
  )

let show_fun _ = "<fun>"
let pp_fun fmt f = pp_of_show show_fun fmt f
let dump_fun fmt f = dump_of_show show_fun fmt f
let display_fun f = display_of_show show_fun f
let print_fun f = print_of_show show_fun f
let prerr_fun f = prerr_of_show show_fun f

let pp_block opening box delim closing fmt pp_elems =
  Format.fprintf fmt "%s" opening;
  Format.fprintf fmt box;
  let open Fstream in
  (match Lazy.force pp_elems with
   | Nil -> ()
   | Cons (pp_elem, pp_elems) ->
     pp_elem fmt;
     Fstream.iter (fun pp -> Format.fprintf fmt "%s@ " delim; pp fmt)
       pp_elems);
  Format.fprintf fmt "@]%s" closing

let pp_seq opening delim closing fmt pp_elems =
  pp_block opening "@[<hov>" delim closing fmt pp_elems

let pp_record pp_fields fmt x =
  let pp_elem (name, f) fmt = Format.fprintf fmt "%s = " name; f fmt x in
  pp_block "{" "@[<hv>" ";" "}" fmt
    (Fstream.map pp_elem (Fstream.of_list pp_fields))

let const k _ = k

let pp_list_ix pp_elem fmt xs =
  pp_seq "[" ";" "]" fmt
    (Fstream.mapi (fun i x fmt -> pp_elem i fmt x) (Fstream.of_list xs))

let pp_list pp_elem fmt xs = pp_list_ix (const pp_elem) fmt xs

(* Darn value restriction.  *)
let dump_list pp_elem = dump_of_pp (pp_list pp_elem)
let show_list pp_elem = show_of_pp (pp_list pp_elem)
let display_list pp_elem = display_of_pp (pp_list pp_elem)
let print_list pp_elem = print_of_pp (pp_list pp_elem)
let prerr_list pp_elem = prerr_of_pp (pp_list pp_elem)

(* Darn value restriction.  *)
let dump_list_ix pp_elem = dump_of_pp (pp_list_ix pp_elem)
let show_list_ix pp_elem = show_of_pp (pp_list_ix pp_elem)
let display_list_ix pp_elem = display_of_pp (pp_list_ix pp_elem)
let print_list_ix pp_elem = print_of_pp (pp_list_ix pp_elem)
let prerr_list_ix pp_elem = prerr_of_pp (pp_list_ix pp_elem)

let pp_array_like length get opening closing pp_elem fmt xs =
  pp_seq opening ";" closing fmt
    (Fstream.map
       (fun i fmt -> pp_elem fmt (get xs i))
       (Fstream.enum 0 (length xs - 1)))

let pp_array pp_elem fmt xs =
  pp_array_like Array.length Array.get "[|" "|]" pp_elem fmt xs

let dump_array pp_elem = dump_of_pp (pp_array pp_elem)
let show_array pp_elem = show_of_pp (pp_array pp_elem)
let display_array pp_elem = display_of_pp (pp_array pp_elem)
let print_array pp_elem = print_of_pp (pp_array pp_elem)
let prerr_array pp_elem = prerr_of_pp (pp_array pp_elem)

let pp_bigarray1 kind layout =
  (* Hmm...this use of Obj.magic is rather grotesque, but without it the kind
     and layout parameters must be both strings, which means the caller can
     confuse which is which.  Obj.magic in this case actually *increases* type
     safety!
     *)
  let kind =
    if Obj.magic kind = Bigarray.float32 then "Bigarray.float32"
    else if Obj.magic kind = Bigarray.float64 then "Bigarray.float64"
    else if Obj.magic kind = Bigarray.complex32 then "Bigarray.complex32"
    else if Obj.magic kind = Bigarray.complex64 then "Bigarray.complex64"
    else if Obj.magic kind = Bigarray.int8_signed then "Bigarray.int8_signed"
    else if Obj.magic kind = Bigarray.int8_unsigned then
      "Bigarray.int8_unsigned"
    else if Obj.magic kind = Bigarray.int16_signed then "Bigarray.int16_signed"
    else if Obj.magic kind = Bigarray.int16_unsigned then
      "Bigarray.int16_unsigned"
    else if Obj.magic kind = Bigarray.int then "Bigarray.int"
    else if Obj.magic kind = Bigarray.int32 then "Bigarray.int32"
    else if Obj.magic kind = Bigarray.int64 then "Bigarray.int64"
    else if Obj.magic kind = Bigarray.nativeint then "Bigarray.nativeint"
    else if Obj.magic kind = Bigarray.char then "Bigarray.char"
    else invalid_arg "Pprint.pp_bigarray1: Unrecognized bigarray kind"
  and layout =
    if Obj.magic layout = Bigarray.c_layout then "Bigarray.c_layout"
    else if Obj.magic layout = Bigarray.fortran_layout then
      "Bigarray.fortran_layout"
    else invalid_arg "Pprint.pp_bigarray1: Unrecognized bigarray layout"
  in
  fun pp_elem fmt xs ->
  if !read_write_invariance
  then pp_array_like Bigarray.Array1.dim Bigarray.Array1.get
        (Printf.sprintf "(Bigarray.Array1.of_array %s %s [|" kind layout)
        "|])" pp_elem fmt xs
  else pp_array_like Bigarray.Array1.dim Bigarray.Array1.get
         "[|" "|]" pp_elem fmt xs

let dump_bigarray1 kind layout pp_elem =
  dump_of_pp (pp_bigarray1 kind layout pp_elem)
let show_bigarray1 kind layout pp_elem =
  show_of_pp (pp_bigarray1 kind layout pp_elem)
let display_bigarray1 kind layout pp_elem =
  display_of_pp (pp_bigarray1 kind layout pp_elem)
let print_bigarray1 kind layout pp_elem =
  print_of_pp (pp_bigarray1 kind layout pp_elem)
let prerr_bigarray1 kind layout pp_elem =
  prerr_of_pp (pp_bigarray1 kind layout pp_elem)

let pp_tuple pp_fields fmt x =
  let pp_elem f fmt = f fmt x in
  pp_seq "(" "," ")" fmt (Fstream.map pp_elem (Fstream.of_list pp_fields))
let show_tuple pp_fields x = show_of_pp (pp_tuple pp_fields) x
let dump_tuple pp_fields = dump_of_pp (pp_tuple pp_fields)
let display_tuple pp_fields = display_of_pp (pp_tuple pp_fields)
let print_tuple pp_fields = print_of_pp (pp_tuple pp_fields)
let prerr_tuple pp_fields = prerr_of_pp (pp_tuple pp_fields)

let pp_pair pp1 pp2 = pp_tuple [(fun fmt (x,y) -> pp1 fmt x);
                                (fun fmt (x,y) -> pp2 fmt y)]
let show_pair pp1 pp2 = show_of_pp (pp_pair pp1 pp2)
let dump_pair pp1 pp2 = dump_of_pp (pp_pair pp1 pp2)
let display_pair pp1 pp2 = display_of_pp (pp_pair pp1 pp2)
let print_pair pp1 pp2 = print_of_pp (pp_pair pp1 pp2)
let prerr_pair pp1 pp2 = prerr_of_pp (pp_pair pp1 pp2)

let pp_triple pp1 pp2 pp3 =
  pp_tuple [(fun fmt (x,y,z) -> pp1 fmt x);
            (fun fmt (x,y,z) -> pp2 fmt y);
            (fun fmt (x,y,z) -> pp3 fmt z)]
let show_triple pp1 pp2 pp3 = show_of_pp (pp_triple pp1 pp2 pp3)
let dump_triple pp1 pp2 pp3 = dump_of_pp (pp_triple pp1 pp2 pp3)
let display_triple pp1 pp2 pp3 = display_of_pp (pp_triple pp1 pp2 pp3)
let print_triple pp1 pp2 pp3 = print_of_pp (pp_triple pp1 pp2 pp3)
let prerr_triple pp1 pp2 pp3 = prerr_of_pp (pp_triple pp1 pp2 pp3)

module PrettyPrim = struct
  let print_string = print_string
  let prerr_string = prerr_string
  let print_float = print_float
  let prerr_float = prerr_float
  let print_char = print_char
  let prerr_char = prerr_char
end
