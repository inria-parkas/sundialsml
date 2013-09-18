(* Utility functions for writing camlp4 extensions.  *)

module Make (Syntax : Camlp4.Sig.Camlp4Syntax) =
struct
  open Syntax
  open Ast

  (* Generic list-processing functions.  *)

  let mapi f xs =
    let rec go acc i = function
      | [] -> List.rev acc
      | x::xs -> go (f i x::acc) (i+1) xs
    in go [] 0 xs

  let map_last f xs =
    match List.rev xs with
    | [] -> invalid_arg "map_last: empty list"
    | x::xs -> List.rev (f x::xs)

  let rec drop_while f = function
    | x::xs when f x -> drop_while f xs
    | xs -> xs

  let rec fold_left_1 f = function
    | [] -> failwith "internal error in camlp4 syntax extension: \
                      fold_left_1 got []"
    | x::xs -> List.fold_left f x xs

  let rec fold_right_1 f xs =
    match List.rev xs with
    | [] -> failwith "internal error in camlp4 syntax extension: \
                      fold_left_1 got []"
    | x::xs -> List.fold_left (fun x y -> f y x) x xs

  (* Printing: idea due to Jeremy Yallop
     https://groups.google.com/forum/#!topic/fa.caml/81RS0IVlS9Y *)

  let string_of_of_format format ?(one_line=true) x =
    let buf = Buffer.create 100 in
    let bfmt = Format.formatter_of_buffer buf in
    format bfmt x;
    Format.pp_print_flush bfmt ();
    let str = Buffer.contents buf in
    (* Collapse to one line, so that it's suitable for an error message.  *)
    if one_line then
      for i = 0 to String.length str - 1 do
        if str.[i] = '\n' then str.[i] <- ' '
      done;
    str

  let print_of_format format x =
    format Format.std_formatter x;
    Format.pp_print_flush Format.std_formatter ()

  let make_printer ?curry_constr ?comments ?semisep () =
    let module P = Camlp4.Printers.OCaml.Make(Syntax) in
    let o = 
      object (self)
        inherit (P.printer ?curry_constr ?comments ())

        method string_of_anti = string_of_of_format self#anti
        method print_anti = print_of_format self#anti
        method string_of_apply_expr = string_of_of_format self#apply_expr
        method print_apply_expr = print_of_format self#apply_expr
        method string_of_binding = string_of_of_format self#binding
        method print_binding = print_of_format self#binding
        method string_of_class_declaration = string_of_of_format self#class_declaration
        method print_class_declaration = print_of_format self#class_declaration
        method string_of_class_expr = string_of_of_format self#class_expr
        method print_class_expr = print_of_format self#class_expr
        method string_of_class_params = string_of_of_format self#class_params
        method print_class_params = print_of_format self#class_params
        method string_of_class_sig_item = string_of_of_format self#class_sig_item
        method print_class_sig_item = print_of_format self#class_sig_item
        method string_of_class_str_item = string_of_of_format self#class_str_item
        method print_class_str_item = print_of_format self#class_str_item
        method string_of_class_type = string_of_of_format self#class_type
        method print_class_type = print_of_format self#class_type
        method string_of_constrain = string_of_of_format self#constrain
        method print_constrain = print_of_format self#constrain
        method string_of_constructor_type = string_of_of_format self#constructor_type
        method print_constructor_type = print_of_format self#constructor_type
        method string_of_ctyp = string_of_of_format self#ctyp
        method print_ctyp = print_of_format self#ctyp
        method string_of_ctyp1 = string_of_of_format self#ctyp1
        method print_ctyp1 = print_of_format self#ctyp1
        method string_of_direction_flag = string_of_of_format self#direction_flag
        method print_direction_flag = print_of_format self#direction_flag
        method string_of_dot_expr = string_of_of_format self#dot_expr
        method print_dot_expr = print_of_format self#dot_expr
        method string_of_expr = string_of_of_format self#expr
        method print_expr = print_of_format self#expr
        method string_of_expr_list = string_of_of_format self#expr_list
        method print_expr_list = print_of_format self#expr_list
        method string_of_fun_binding = string_of_of_format self#fun_binding
        method print_fun_binding = print_of_format self#fun_binding
        method string_of_functor_arg = string_of_of_format self#functor_arg
        method print_functor_arg = print_of_format self#functor_arg
        method string_of_functor_args = string_of_of_format self#functor_args
        method print_functor_args = print_of_format self#functor_args
        method string_of_ident = string_of_of_format self#ident
        method print_ident = print_of_format self#ident
        method string_of_implem = string_of_of_format self#implem
        method print_implem = print_of_format self#implem
        method string_of_interf = string_of_of_format self#interf
        method print_interf = print_of_format self#interf
        method string_of_match_case = string_of_of_format self#match_case
        method print_match_case = print_of_format self#match_case
        method string_of_match_case_aux = string_of_of_format self#match_case_aux
        method print_match_case_aux = print_of_format self#match_case_aux
        method string_of_module_expr = string_of_of_format self#module_expr
        method print_module_expr = print_of_format self#module_expr
        method string_of_module_rec_binding = string_of_of_format self#module_rec_binding
        method print_module_rec_binding = print_of_format self#module_rec_binding
        method string_of_module_type = string_of_of_format self#module_type
        method print_module_type = print_of_format self#module_type
        method string_of_mutable_flag = string_of_of_format self#mutable_flag
        method print_mutable_flag = print_of_format self#mutable_flag
        method string_of_numeric = string_of_of_format self#numeric
        method print_numeric = print_of_format self#numeric
        method string_of_override_flag = string_of_of_format self#override_flag
        method print_override_flag = print_of_format self#override_flag
        method string_of_patt = string_of_of_format self#patt
        method print_patt = print_of_format self#patt
        method string_of_patt1 = string_of_of_format self#patt1
        method print_patt1 = print_of_format self#patt1
        method string_of_patt2 = string_of_of_format self#patt2
        method print_patt2 = print_of_format self#patt2
        method string_of_patt3 = string_of_of_format self#patt3
        method print_patt3 = print_of_format self#patt3
        method string_of_patt4 = string_of_of_format self#patt4
        method print_patt4 = print_of_format self#patt4
        method string_of_patt5 = string_of_of_format self#patt5
        method print_patt5 = print_of_format self#patt5
        method string_of_patt_class_expr_fun_args = string_of_of_format self#patt_class_expr_fun_args
        method print_patt_class_expr_fun_args = print_of_format self#patt_class_expr_fun_args
        method string_of_patt_expr_fun_args = string_of_of_format self#patt_expr_fun_args
        method print_patt_expr_fun_args = print_of_format self#patt_expr_fun_args
        method string_of_patt_tycon = string_of_of_format self#patt_tycon
        method print_patt_tycon = print_of_format self#patt_tycon
        method string_of_private_flag = string_of_of_format self#private_flag
        method print_private_flag = print_of_format self#private_flag
        method string_of_quoted_string = string_of_of_format self#quoted_string
        method print_quoted_string = print_of_format self#quoted_string
        method string_of_raise_match_failure = string_of_of_format self#raise_match_failure
        method print_raise_match_failure = print_of_format self#raise_match_failure
        method string_of_rec_flag = string_of_of_format self#rec_flag
        method print_rec_flag = print_of_format self#rec_flag
        method string_of_record_binding = string_of_of_format self#record_binding
        method print_record_binding = print_of_format self#record_binding
        method string_of_seq = string_of_of_format self#seq
        method print_seq = print_of_format self#seq
        method string_of_sig_item = string_of_of_format self#sig_item
        method print_sig_item = print_of_format self#sig_item
        method string_of_simple_ctyp = string_of_of_format self#simple_ctyp
        method print_simple_ctyp = print_of_format self#simple_ctyp
        method string_of_simple_expr = string_of_of_format self#simple_expr
        method print_simple_expr = print_of_format self#simple_expr
        method string_of_simple_module_expr = string_of_of_format self#simple_module_expr
        method print_simple_module_expr = print_of_format self#simple_module_expr
        method string_of_simple_patt = string_of_of_format self#simple_patt
        method print_simple_patt = print_of_format self#simple_patt
        method string_of_str_item = string_of_of_format self#str_item
        method print_str_item = print_of_format self#str_item
        method string_of_string = string_of_of_format self#string
        method print_string = print_of_format self#string
        method string_of_sum_type = string_of_of_format self#sum_type
        method print_sum_type = print_of_format self#sum_type
        method string_of_type_params = string_of_of_format self#type_params
        method print_type_params = print_of_format self#type_params
        method string_of_var = string_of_of_format self#var
        method print_var = print_of_format self#var
        method string_of_virtual_flag = string_of_of_format self#virtual_flag
        method print_virtual_flag = print_of_format self#virtual_flag
        method string_of_with_constraint = string_of_of_format self#with_constraint
        method print_with_constraint = print_of_format self#with_constraint
      end
    in
    match semisep with
    | None -> o
    | Some sep -> o#set_semisep sep

  let ocaml_printer = make_printer ()

  let string_of_anti = string_of_of_format ocaml_printer#anti
  let print_anti = print_of_format ocaml_printer#anti
  let string_of_apply_expr = string_of_of_format ocaml_printer#apply_expr
  let print_apply_expr = print_of_format ocaml_printer#apply_expr
  let string_of_binding = string_of_of_format ocaml_printer#binding
  let print_binding = print_of_format ocaml_printer#binding
  let string_of_class_declaration = string_of_of_format ocaml_printer#class_declaration
  let print_class_declaration = print_of_format ocaml_printer#class_declaration
  let string_of_class_expr = string_of_of_format ocaml_printer#class_expr
  let print_class_expr = print_of_format ocaml_printer#class_expr
  let string_of_class_params = string_of_of_format ocaml_printer#class_params
  let print_class_params = print_of_format ocaml_printer#class_params
  let string_of_class_sig_item = string_of_of_format ocaml_printer#class_sig_item
  let print_class_sig_item = print_of_format ocaml_printer#class_sig_item
  let string_of_class_str_item = string_of_of_format ocaml_printer#class_str_item
  let print_class_str_item = print_of_format ocaml_printer#class_str_item
  let string_of_class_type = string_of_of_format ocaml_printer#class_type
  let print_class_type = print_of_format ocaml_printer#class_type
  let string_of_constrain = string_of_of_format ocaml_printer#constrain
  let print_constrain = print_of_format ocaml_printer#constrain
  let string_of_constructor_type = string_of_of_format ocaml_printer#constructor_type
  let print_constructor_type = print_of_format ocaml_printer#constructor_type
  let string_of_ctyp = string_of_of_format ocaml_printer#ctyp
  let print_ctyp = print_of_format ocaml_printer#ctyp
  let string_of_ctyp1 = string_of_of_format ocaml_printer#ctyp1
  let print_ctyp1 = print_of_format ocaml_printer#ctyp1
  let string_of_direction_flag = string_of_of_format ocaml_printer#direction_flag
  let print_direction_flag = print_of_format ocaml_printer#direction_flag
  let string_of_dot_expr = string_of_of_format ocaml_printer#dot_expr
  let print_dot_expr = print_of_format ocaml_printer#dot_expr
  let string_of_expr = string_of_of_format ocaml_printer#expr
  let print_expr = print_of_format ocaml_printer#expr
  let string_of_expr_list = string_of_of_format ocaml_printer#expr_list
  let print_expr_list = print_of_format ocaml_printer#expr_list
  let string_of_fun_binding = string_of_of_format ocaml_printer#fun_binding
  let print_fun_binding = print_of_format ocaml_printer#fun_binding
  let string_of_functor_arg = string_of_of_format ocaml_printer#functor_arg
  let print_functor_arg = print_of_format ocaml_printer#functor_arg
  let string_of_functor_args = string_of_of_format ocaml_printer#functor_args
  let print_functor_args = print_of_format ocaml_printer#functor_args
  let string_of_ident = string_of_of_format ocaml_printer#ident
  let print_ident = print_of_format ocaml_printer#ident
  let string_of_implem = string_of_of_format ocaml_printer#implem
  let print_implem = print_of_format ocaml_printer#implem
  let string_of_interf = string_of_of_format ocaml_printer#interf
  let print_interf = print_of_format ocaml_printer#interf
  let string_of_match_case = string_of_of_format ocaml_printer#match_case
  let print_match_case = print_of_format ocaml_printer#match_case
  let string_of_match_case_aux = string_of_of_format ocaml_printer#match_case_aux
  let print_match_case_aux = print_of_format ocaml_printer#match_case_aux
  let string_of_module_expr = string_of_of_format ocaml_printer#module_expr
  let print_module_expr = print_of_format ocaml_printer#module_expr
  let string_of_module_rec_binding = string_of_of_format ocaml_printer#module_rec_binding
  let print_module_rec_binding = print_of_format ocaml_printer#module_rec_binding
  let string_of_module_type = string_of_of_format ocaml_printer#module_type
  let print_module_type = print_of_format ocaml_printer#module_type
  let string_of_mutable_flag = string_of_of_format ocaml_printer#mutable_flag
  let print_mutable_flag = print_of_format ocaml_printer#mutable_flag
  let string_of_numeric = string_of_of_format ocaml_printer#numeric
  let print_numeric = print_of_format ocaml_printer#numeric
  let string_of_override_flag = string_of_of_format ocaml_printer#override_flag
  let print_override_flag = print_of_format ocaml_printer#override_flag
  let string_of_patt = string_of_of_format ocaml_printer#patt
  let print_patt = print_of_format ocaml_printer#patt
  let string_of_patt1 = string_of_of_format ocaml_printer#patt1
  let print_patt1 = print_of_format ocaml_printer#patt1
  let string_of_patt2 = string_of_of_format ocaml_printer#patt2
  let print_patt2 = print_of_format ocaml_printer#patt2
  let string_of_patt3 = string_of_of_format ocaml_printer#patt3
  let print_patt3 = print_of_format ocaml_printer#patt3
  let string_of_patt4 = string_of_of_format ocaml_printer#patt4
  let print_patt4 = print_of_format ocaml_printer#patt4
  let string_of_patt5 = string_of_of_format ocaml_printer#patt5
  let print_patt5 = print_of_format ocaml_printer#patt5
  let string_of_patt_class_expr_fun_args = string_of_of_format ocaml_printer#patt_class_expr_fun_args
  let print_patt_class_expr_fun_args = print_of_format ocaml_printer#patt_class_expr_fun_args
  let string_of_patt_expr_fun_args = string_of_of_format ocaml_printer#patt_expr_fun_args
  let print_patt_expr_fun_args = print_of_format ocaml_printer#patt_expr_fun_args
  let string_of_patt_tycon = string_of_of_format ocaml_printer#patt_tycon
  let print_patt_tycon = print_of_format ocaml_printer#patt_tycon
  let string_of_private_flag = string_of_of_format ocaml_printer#private_flag
  let print_private_flag = print_of_format ocaml_printer#private_flag
  let string_of_quoted_string = string_of_of_format ocaml_printer#quoted_string
  let print_quoted_string = print_of_format ocaml_printer#quoted_string
  let string_of_raise_match_failure = string_of_of_format ocaml_printer#raise_match_failure
  let print_raise_match_failure = print_of_format ocaml_printer#raise_match_failure
  let string_of_rec_flag = string_of_of_format ocaml_printer#rec_flag
  let print_rec_flag = print_of_format ocaml_printer#rec_flag
  let string_of_record_binding = string_of_of_format ocaml_printer#record_binding
  let print_record_binding = print_of_format ocaml_printer#record_binding
  let string_of_seq = string_of_of_format ocaml_printer#seq
  let print_seq = print_of_format ocaml_printer#seq
  let string_of_sig_item = string_of_of_format ocaml_printer#sig_item
  let print_sig_item = print_of_format ocaml_printer#sig_item
  let string_of_simple_ctyp = string_of_of_format ocaml_printer#simple_ctyp
  let print_simple_ctyp = print_of_format ocaml_printer#simple_ctyp
  let string_of_simple_expr = string_of_of_format ocaml_printer#simple_expr
  let print_simple_expr = print_of_format ocaml_printer#simple_expr
  let string_of_simple_module_expr = string_of_of_format ocaml_printer#simple_module_expr
  let print_simple_module_expr = print_of_format ocaml_printer#simple_module_expr
  let string_of_simple_patt = string_of_of_format ocaml_printer#simple_patt
  let print_simple_patt = print_of_format ocaml_printer#simple_patt
  let string_of_str_item = string_of_of_format ocaml_printer#str_item
  let print_str_item = print_of_format ocaml_printer#str_item
  let string_of_string = string_of_of_format ocaml_printer#string
  let print_string = print_of_format ocaml_printer#string
  let string_of_sum_type = string_of_of_format ocaml_printer#sum_type
  let print_sum_type = print_of_format ocaml_printer#sum_type
  let string_of_type_params = string_of_of_format ocaml_printer#type_params
  let print_type_params = print_of_format ocaml_printer#type_params
  let string_of_var = string_of_of_format ocaml_printer#var
  let print_var = print_of_format ocaml_printer#var
  let string_of_virtual_flag = string_of_of_format ocaml_printer#virtual_flag
  let print_virtual_flag = print_of_format ocaml_printer#virtual_flag
  let string_of_with_constraint = string_of_of_format ocaml_printer#with_constraint
  let print_with_constraint = print_of_format ocaml_printer#with_constraint

  let hashtbl_of_list xs =
    let h = Hashtbl.create 10 in
    ignore (List.fold_left (fun i x -> Hashtbl.add h x i; i+1) 0 xs);
    h

  let hashtbl_of_assocs xs =
    let h = Hashtbl.create 10 in
    List.iter (fun (x,y) -> Hashtbl.add h x y) xs;
    h

  let fail_at _loc msg = Loc.raise _loc (Failure msg)
  let fail_ctyp t msg = fail_at (loc_of_ctyp t) (msg ^ ": " ^ string_of_ctyp t)

  let string_list_of_ident id =
    let go = function
      | <:ident<$lid:x$>> -> x
      | <:ident<$uid:x$>> -> x
      | id -> fail_at (loc_of_ident id) ("unrecognized component in \
                                          identifier: " ^ string_of_ident id)
    in List.map go (list_of_ident id [])

  (* Smart term constructors *)

  let smart_comb ctor e e' =
    match e, e' with
    | ExNil _loc, ExNil _loc' -> ExNil (Loc.merge _loc _loc')
    | e, ExNil _ -> e
    | ExNil _, e' -> e'
    | e, e' -> ctor (Loc.merge (loc_of_expr e) (loc_of_expr e')) e e'

  let sem = smart_comb (fun l e e' -> ExSem (l, e, e'))
  let com = smart_comb (fun l e e' -> ExCom (l, e, e'))
  let app = smart_comb (fun l e e' -> ExApp (l, e, e'))

  let smart_comb_patt ctor e e' =
    match e, e' with
    | PaNil _loc, PaNil _loc' -> PaNil (Loc.merge _loc _loc')
    | e, PaNil _ -> e
    | PaNil _, e' -> e'
    | e, e' -> ctor (Loc.merge (loc_of_patt e) (loc_of_patt e')) e e'

  let sem_patt = smart_comb_patt (fun l e e' -> PaSem (l, e, e'))
  let com_patt = smart_comb_patt (fun l e e' -> PaCom (l, e, e'))

  let smart_comb_ctyp ctor e e' =
    match e, e' with
    | TyNil _loc, TyNil _loc' -> TyNil (Loc.merge _loc _loc')
    | e, TyNil _ -> e
    | TyNil _, e' -> e'
    | e, e' -> ctor (Loc.merge (loc_of_ctyp e) (loc_of_ctyp e')) e e'

  let sem_ctyp = smart_comb_ctyp (fun l e e' -> TySem (l, e, e'))
  let com_ctyp = smart_comb_ctyp (fun l e e' -> TyCom (l, e, e'))

  let smart_comb_binding ctor e e' =
    match e, e' with
    | BiNil _loc, BiNil _loc' -> BiNil (Loc.merge _loc _loc')
    | e, BiNil _ -> e
    | BiNil _, e' -> e'
    | e, e' -> ctor (Loc.merge (loc_of_binding e) (loc_of_binding e')) e e'

  let and_binding = smart_comb_binding (fun l e e' -> BiAnd (l, e, e'))

  let smart_comb_rec_binding ctor e e' =
    match e, e' with
    | RbNil _loc, RbNil _loc' -> RbNil (Loc.merge _loc _loc')
    | e, RbNil _ -> e
    | RbNil _, e' -> e'
    | e, e' -> ctor (Loc.merge (loc_of_rec_binding e) (loc_of_rec_binding e'))
                 e e'

  let sem_rec_binding = smart_comb_rec_binding (fun l e e' -> RbSem (l, e, e'))

  let id_of_string_list _loc strs =
    let f s =
      if String.length s <= 0 then
        invalid_arg "id_of_string_list found empty string";
      if s.[0] = Char.uppercase s.[0]
      then <:ident<$uid:s$>>
      else <:ident<$lid:s$>>
    in
    idAcc_of_list (List.map f strs)

  (* Mass-term constructors *)

  let fold_expr_with ctor _loc e es = List.fold_left (smart_comb ctor) e es

  let fold_right_expr_with ctor _loc es e =
    List.fold_right (smart_comb ctor) es e

  let fold_patt_with ctor _loc p ps =
    let go p p' =
      match p, p' with
      | PaNil _loc, PaNil _loc' -> PaNil (Loc.merge _loc _loc')
      | p, PaNil _ -> p
      | PaNil _, p' -> p'
      | p, p' -> ctor (Loc.merge (loc_of_patt p) (loc_of_patt p')) p p'
    in List.fold_left go p ps

  let fold_ctyp_with ctor _loc t ts =
    let go t t' =
      match t, t' with
      | TyNil _loc, TyNil _loc' -> TyNil (Loc.merge _loc _loc')
      | t, TyNil _ -> t
      | TyNil _, t' -> t'
      | t, t' -> ctor (Loc.merge (loc_of_ctyp t) (loc_of_ctyp t')) t t'
    in List.fold_left go t ts

  let fold_right_ctyp_with ctor _loc ts t =
    let go t t' =
      match t, t' with
      | TyNil _loc, TyNil _loc' -> TyNil (Loc.merge _loc _loc')
      | t, TyNil _ -> t
      | TyNil _, t' -> t'
      | t, t' -> ctor (Loc.merge (loc_of_ctyp t) (loc_of_ctyp t')) t t'
    in List.fold_right go ts t

  let fold_app f xs =
    match f with
    | ExNil _loc -> fail_at _loc "internal error in deriving: fold_app got \
                                  ExNil as function"
    | _ -> fold_expr_with (fun _loc e e' -> ExApp (_loc, e, e'))
             (loc_of_expr f) f xs

  let fold_app_ctyp f args =
    fold_ctyp_with (fun _loc a t -> TyApp (_loc, a, t))
      (loc_of_ctyp f) f args

  let fold_fun _loc args body =
    let rec go body = function
      | [] -> body
      | arg::args -> go <:expr<fun $lid:arg$ -> $body$>> args
    in go body (List.rev args)

  let fold_arr args body =
    fold_right_ctyp_with (fun _loc t u -> <:ctyp<$t$ -> $u$>>)
      (loc_of_ctyp body) args body

  let forall _loc args body =
    match List.map (fun x -> <:ctyp<'$lid:x$>>) args with
    | [] -> body
    | [arg] -> TyPol (_loc, arg, body)
    | arg::args -> TyPol (_loc, fold_app_ctyp arg args, body)

  let fold_sem _loc es =
    match es with
    | [] -> ExNil _loc
    | e::es -> List.fold_left sem e es

  let fold_com _loc es =
    match es with
    | [] -> ExNil _loc
    | e::es ->
      let go e e' =
        match e, e' with
        | ExNil _loc, ExNil _loc' -> ExNil (Loc.merge _loc _loc')
        | e, ExNil _ -> e
        | ExNil _, e' -> e'
        | e, e' -> ExCom (Loc.merge (loc_of_expr e) (loc_of_expr e'), e, e')
      in List.fold_left go e es

  (* [p1;p2;...] -> <:patt<$p1$, $p2$, ...>> *)
  let fold_app_patt p ps =
    fold_patt_with (fun _loc p p' -> PaApp (_loc, p, p')) (loc_of_patt p) p ps

  (* [e] -> e
     [e1;e2;...] -> <:expr<($e1$, $e2$, ...)>> *)
  let fold_tuple = function
    | [] -> failwith  "error in camlp4 syntax manipulation: fold_tuple []"
    | [e] -> e
    | e::es -> let e = List.fold_left com e es in
               ExTup (loc_of_expr e, e)

  let fold_tuple_patt = function
    | [] -> failwith "error in camlp4 syntax manipulation: fold_tuple_patt []"
    | [p] -> p
    | p::ps -> let p = List.fold_left com_patt p ps in
               PaTup (loc_of_patt p, p)

  (* [p1;p2;...] -> <:binding<$p1$ and $p2$ and ...>> *)
  let fold_and_binding bs =
    match bs with
    | [] -> failwith "error in camlp4 syntax manipulation: \
                      fold_sem_rec_binding []"
    | b::bs -> List.fold_left and_binding b bs

  let fold_sem_rec_binding = function
    | [] -> failwith "error in camlp4 syntax manipulation: \
                      fold_sem_rec_binding []"
    | b::bs -> List.fold_left sem_rec_binding b bs

  (* [m1;m2;...] -> <:match_case<$m1$ | $m2$ | ...>> *)
  let fold_or_match_case _loc ms =
    match ms with
    | [] -> McNil _loc
    | m::ms ->
      let go m m' =
        let _loc = Loc.merge (loc_of_match_case m) (loc_of_match_case m') in
        <:match_case<$m$ | $m'$>>
      in List.fold_left go m ms

  let fold_str_item _loc ss =
    let go s s' =
      match s, s' with
      | StNil _loc, StNil _loc' -> StNil (Loc.merge _loc _loc')
      | s, StNil _ -> s
      | StNil _, s' -> s'
      | s, s' -> StSem (Loc.merge (loc_of_str_item s) (loc_of_str_item s'),
                        s, s')
    in List.fold_right go ss (StNil _loc)

  (* Auxiliary functions for reducing the nesting depth of antiquotes.  *)

  let meta_lid _loc name = <:expr<$lid:name$>>
  let meta_lid_ident _loc name = <:ident<$lid:name$>>
  let meta_lid_patt _loc name = <:patt<$lid:name$>>

  let meta_ident ident =
    let _loc = loc_of_ident ident in
    <:expr<$id:ident$>>

  let meta_uid _loc name = <:expr<$uid:name$>>
  let meta_uid_ident _loc name = <:ident<$uid:name$>>
  let meta_uid_patt _loc name = <:patt<$uid:name$>>

  (** [meta_com : <<'a>> -> <<'b>> -> <<'a * 'b>>]
      [meta_com e1 e2 = <:expr< <:expr<$$e1$$, $$e2$$>> >>]
      where $$ is to be interpreted as nested antiquotes and where
      [<'a>] denotes expr of type ['a] and [*] means [ExCom], not [ExTup].

      For example, if [f] is a function at stage-1 (i.e. one that exists in the
      generated code) of type [int -> expr], then
      [meta_com <:expr< f 0 >> <:expr< f 1 >> =
         <:expr< <:expr<f 0, f 1>> >>]
  *)
  let meta_com locvar e1 e2 =
    let _loc = Loc.merge (loc_of_expr e1) (loc_of_expr e2) in
    <:expr< Ast.ExCom ($locvar$, $e1$, $e2$) >>

  let meta_fold_com locvar = function
    | [] -> failwith "error in camlp4 syntax manipulation: meta_fold_com []"
    | e::es -> List.fold_left (meta_com locvar) e es

  (* [e1;e2;...] -> <:expr< <:expr<($$e1$$, $$e2$$, ...)>> >> *)
  let meta_tuple _loc locvar = function
    | [] -> <:expr< <:expr<()>> >>
    | es -> <:expr< Ast.ExTup ($locvar$, $meta_fold_com locvar es$)>>

  let meta_app locvar e1 e2 =
    let _loc = Loc.merge (loc_of_expr e1) (loc_of_expr e2) in
    <:expr< Ast.ExApp ($locvar$, $e1$, $e2$) >>

  let meta_fold_app locvar = function
    | [] -> failwith "error in camlp4 syntax manipulation: meta_fold_com []"
    | e::es -> List.fold_left (meta_app locvar) e es

  let meta_sem_rec_binding locvar e1 e2 =
    let _loc = Loc.merge (loc_of_expr e1) (loc_of_expr e2) in
    <:expr< Ast.RbSem ($locvar$, $e1$, $e2$) >>

  let meta_fold_sem_rec_binding locvar = function
    | [] -> failwith "error in camlp4 syntax manipulation: meta_fold_com []"
    | e::es -> List.fold_left (meta_sem_rec_binding locvar) e es

end
