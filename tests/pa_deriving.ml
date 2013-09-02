(* A camlp4 extension for generating pretty printers for data types.  Requires
OCaml 3.12 or newer.  The code

  open Pprint
  type 'a foo = Foo of ('a, int) bar
  deriving pretty

generates

  type 'a foo = Foo of ('a, int) bar

  let pp_foo pp0 ?(prec=0) fmt = function
    | Foo x -> pp_string_noquote fmt "Foo ";

Note that open Pprint is strictly necessary.  If the definition of foo refers
to another type baz, for example, this generator will call a pretty-printer
named pp_baz.  The user can either write pp_baz by hand or derive it with
another deriving clause at the definition of type baz, but in any case the user
is responsible for ensuring that pp_baz exists and is visible.  Pprint supplies
these functions for base types like int, bool, etc., and it also supplies some
auxiliary functions that this generator relies on.

*)

module Id : Camlp4.Sig.Id = struct
  let name = "pa_deriving_pprint"
  let version = "1.0"
end

module Make (Syntax : Camlp4.Sig.Camlp4Syntax) =
struct
  open Camlp4.Sig
  include Syntax
  open Ast
  (*open Camlp4aux
  module GenHelper = GenHelper (Syntax)*)

  (*
   * Auxiliary defs
   *)

  (* Technique due to Jeremy Yallop
     https://groups.google.com/forum/#!topic/fa.caml/81RS0IVlS9Y *)
  let ocaml_printer =
    let module P = Camlp4.Printers.OCaml.Make(Syntax) in
    new P.printer ()

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

  let string_of_str_item = string_of_of_format ocaml_printer#str_item
  let string_of_ctyp = string_of_of_format ocaml_printer#ctyp
  let string_of_expr = string_of_of_format ocaml_printer#expr
  let string_of_patt = string_of_of_format ocaml_printer#patt
  let string_of_ident = string_of_of_format ocaml_printer#ident
  let string_of_binding = string_of_of_format ocaml_printer#binding

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

  let mapi f xs =
    let rec go acc i = function
      | [] -> List.rev acc
      | x::xs -> go (f i x::acc) (i+1) xs
    in go [] 0 xs

  let rec drop_while f = function
    | x::xs when f x -> drop_while f xs
    | xs -> xs

  let string_list_of_ident id =
    let go = function
      | <:ident<$lid:x$>> -> x
      | <:ident<$uid:x$>> -> x
      | id -> fail_at (loc_of_ident id) ("unrecognized component in \
                                          identifier: " ^ string_of_ident id)
    in List.map go (list_of_ident id [])

  (* Mass-term constructors *)

  let fold_expr_with ctor _loc e es =
    let go e e' =
      match e, e' with
      | ExNil _loc, ExNil _loc' -> ExNil (Loc.merge _loc _loc')
      | e, ExNil _ -> e
      | ExNil _, e' -> e'
      | e, e' -> ctor (Loc.merge (loc_of_expr e) (loc_of_expr e')) e e'
    in List.fold_left go e es

  let fold_right_expr_with ctor _loc es e =
    let go e e' =
      match e, e' with
      | ExNil _loc, ExNil _loc' -> ExNil (Loc.merge _loc _loc')
      | e, ExNil _ -> e
      | ExNil _, e' -> e'
      | e, e' -> ctor (Loc.merge (loc_of_expr e) (loc_of_expr e')) e e'
    in List.fold_right go es e

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

  let fold_semi _loc es =
    match es with
    | [] -> ExNil _loc
    | e::es ->
      let go e e' =
        match e, e' with
        | ExNil _loc, ExNil _loc' -> ExNil (Loc.merge _loc _loc')
        | e, ExNil _ -> e
        | ExNil _, e' -> e'
        | e, e' -> ExSem (Loc.merge (loc_of_expr e) (loc_of_expr e'), e, e')
      in List.fold_left go e es

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
  let fold_com_patt _loc ps =
    match ps with
    | [] -> PaNil _loc
    | p::ps -> fold_patt_with (fun _loc p p' -> PaCom (_loc, p, p')) _loc p ps

  (* [p1;p2;...] -> <:patt<$p1$, $p2$, ...>> *)
  let fold_app_patt p ps =
    fold_patt_with (fun _loc p p' -> PaApp (_loc, p, p')) (loc_of_patt p) p ps

  let fold_tuple_patt = function
    | [] | [_] -> failwith "internal error in deriving: fold_tuple_patt \
                            argument list too short"
    | ps -> let p = fold_com_patt Loc.ghost ps in
            PaTup (loc_of_patt p, p)

  (* [p1;p2;...] -> <:binding<$p1$ and $p2$ and ...>> *)
  let fold_and_binding _loc bs =
    match bs with
    | [] -> BiNil _loc
    | b::bs ->
      let go b b' =
        let _loc = Loc.merge (loc_of_binding b) (loc_of_binding b') in
        <:binding<$b$ and $b'$>>
      in List.fold_left go b bs

  (* [m1;m2;...] -> <:match_case<$m1$ | $m2$ | ...>> *)
  let fold_or_match_case _loc ms =
    match ms with
    | [] -> McNil _loc
    | m::ms ->
      let go m m' =
        let _loc = Loc.merge (loc_of_match_case m) (loc_of_match_case m') in
        <:match_case<$m$ | $m'$>>
      in List.fold_left go m ms

  (* Simplified representation of a type declaration.  *)
  type type_decl_simpl =
    (* list of (loc, type parameters, type constructor, rhs)
       Type parameters are listed in the order they appear in the source code
       which is opposite from the order in which they appear in the AST.  *)
      TDecl of Loc.t * (Loc.t * string list * string * type_def_simpl) list
  and  type_def_simpl =
    (* list of (field name, mutability, field type) *)
    | TRecord of Loc.t * (Loc.t * string * bool * type_expr_simpl) list
    (* list of (constructor name, list of arguments) *)
    | TVariant of Loc.t * (Loc.t * string * type_expr_simpl list) list
    | TAlias of type_expr_simpl
  and  type_expr_simpl =
    | TApp of Loc.t * ident * type_expr_simpl list (* fun, args (may be []) *)
    | TTup of Loc.t * type_expr_simpl list
    | TFun of Loc.t * type_expr_simpl * type_expr_simpl
    | TVar of Loc.t * string

  let rec type_decl_simpl_of_ctyp ctyp =
    let rec go ctyp eqs =
      match ctyp with
      | <:ctyp@_loc<$d1$ and $d2$>> -> go d1 (go d2 eqs)
      | TyDcl (_loc, tctor, params, rhs, []) ->
        (_loc, List.map tyquo_name params, tctor, type_def_simpl_of_ctyp rhs)
        ::eqs
      | t -> fail_ctyp t "deriving can't handle this type declaration"
    in TDecl (loc_of_ctyp ctyp, go ctyp [])

  (* (AST of 'foo) -> "foo" *)
  and tyquo_name = function
    | TyQuo (_, name) -> name
    | t -> fail_ctyp t "deriving pretty can't handle this type parameter"

  and type_def_simpl_of_ctyp = function
    | <:ctyp@_loc<| $t$ >> -> TVariant (_loc, flatten_variants t)
    | <:ctyp@_loc<{ $t$ }>> -> TRecord (_loc, flatten_record t)
    | t -> TAlias (type_expr_simpl_of_ctyp t)

  and flatten_variants t =
    let rec go t cases =
      match t with
      | <:ctyp<$t$ | $u$>> -> go t (go u cases)
      | <:ctyp@_loc<$uid:ctor$>> -> (_loc, ctor, [])::cases
      | <:ctyp@_loc<$uid:ctor$ of $t$>> -> (_loc, ctor, go_of t [])::cases
      | t -> fail_ctyp t "deriving: I didn't expect this in a variant type \
                          declaration, I'm confused"
    and go_of t args =
      match t with
      | <:ctyp<$t$ and $u$>> -> go_of t (go_of u args)
      | t -> type_expr_simpl_of_ctyp t::args
    in go t []

  and flatten_record t =
    let rec go t fields =
      match t with
      | TyCol (_loc, <:ctyp<$lid:f$>>, <:ctyp<mutable $t$>>) ->
        (_loc, f, true, type_expr_simpl_of_ctyp t)::fields
      | TyCol (_loc, <:ctyp<$lid:f$>>, <:ctyp<$t$>>) ->
        (_loc, f, false, type_expr_simpl_of_ctyp t)::fields
      | <:ctyp<$f$; $g$>> -> go f (go g fields)
      | t -> fail_ctyp t "deriving: I didn't expect this in a record's \
                          field declaration, I'm confused"
    in go t []

  and type_expr_simpl_of_ctyp = function
    | <:ctyp@_loc<'$lid:x$>> -> TVar (_loc, x)
    | <:ctyp@_loc<$id:x$>> -> TApp (_loc, x, [])
    | <:ctyp@_loc<$_$ $_$>> as t -> let (f, args) = flatten_tapp t in
                                    TApp (_loc, f, args)
    | <:ctyp@_loc<$t$ -> $u$>> -> TFun (_loc,
                                        type_expr_simpl_of_ctyp t,
                                        type_expr_simpl_of_ctyp u)
    | TyTup (_loc, t) -> TTup (_loc, flatten_tuple t)
    | t -> fail_ctyp t "deriving: can't handle this kind of type expression"

  and flatten_tapp t0 =
    (* NB: (a,b,c) foo is parsed as (c (b (a foo))).  *)
    let rec go t args =
      match t with
      | <:ctyp<$arg$ $f$>> -> go f (type_expr_simpl_of_ctyp arg::args)
      | <:ctyp<$id:x$>> -> (x, args)
      | t -> fail_ctyp t0 "deriving: I don't know what to do with this type"
    in go t0 []

  and flatten_tuple t =
    let rec go t components =
      match t with
      | TySta (_loc, t, u) -> go t (go u components)
      | t -> type_expr_simpl_of_ctyp t::components
    in go t []

  let rec flatten_longident = function
    | <:ident<$lid:x$>> -> x
    | <:ident<$uid:x$>> -> x
    | <:ident<$a$.$b$>> -> flatten_longident a ^ "_" ^ flatten_longident b
    | id -> fail_at (loc_of_ident id) ("deriving can't handle this kind of \
                                        identifier")

  (* Pretty-printer derivation *)

  (* Generate a type for the pretty-printer of the given type.  This is
     necessary since the printer for e.g.
       type 'a foo = ...
       and bar = int foo
     requires polymorphic recursion.  This is the part that needs
     OCaml >= 3.12.  *)
  let pretty_type _loc param_names tctor =
    let params = List.map (fun x -> <:ctyp<'$lid:x$>>) param_names in
    forall _loc param_names
      (fold_arr (List.map (fun x -> <:ctyp<$x$ pp>>) params)
         (TyApp (_loc, <:ctyp<pp>>,
                 fold_app_ctyp <:ctyp<$lid:tctor$>> params)))

  (* Generate a pretty-printer's body.  The resulting expression has type
     'a pp.  *)
  let pretty_body _loc aliases long_prefix omit params rhs =
    (* Table from type parameter name to pp0, pp1, etc. *)
    let pp_table = Hashtbl.create 10 in
    List.iter (fun x -> Hashtbl.add pp_table x
                  <:expr<$lid:"pp"^string_of_int (Hashtbl.length pp_table)$>>)
      params;
    let short_prefix = drop_while (Hashtbl.mem omit) long_prefix in
    let qualify str =
      let prepend mods s = String.concat "." (mods @ [s]) in
      if short_prefix = long_prefix then
        <:expr<$`str:prepend short_prefix str$>>
      else
        <:expr<if !read_write_invariance
               then $`str:prepend long_prefix str$
               else $`str:prepend short_prefix str$>>
    in
    let rec go_type_expr = function
      | TApp (_loc, <:ident<$x$>>, args) ->
        let x = string_list_of_ident x in
        let x = try Hashtbl.find aliases x with Not_found -> x in
        let x = String.concat "_" x in
        fold_app <:expr<$lid:"pp_"^x$>> (List.map go_type_expr args)
      | TTup (_loc, []) -> fail_at _loc
                           "internal error in deriving: TTup (_,[])"
      | TTup (_loc, [x]) -> fail_at _loc
                           "internal error in deriving: TTup (_,[_])"
      | TTup (_loc, [x;y]) -> <:expr<pp_pair $go_type_expr x$
                                             $go_type_expr y$>>
      | TTup (_loc, [x;y;z]) -> <:expr<pp_triple $go_type_expr x$
                                                 $go_type_expr y$
                                                 $go_type_expr z$>>
      | TTup (_loc, xs) ->
        let mk_pat i _ = <:patt<$lid:"x"^string_of_int i$>>
        and mk_pp i x = <:expr<pp_fixarg $go_type_expr x$
                                         $lid:"x"^string_of_int i$>>
        in
        <:expr<fun ?(prec=0) fmt -> function
          | $fold_tuple_patt (mapi mk_pat xs)$ ->
            pp_seq "(@[<hv>" "," "@])" [| $fold_semi _loc (mapi mk_pp xs)$ |]
        >>
      | TVar (_loc, a) ->
        (try Hashtbl.find pp_table a
         with Not_found -> fail_at _loc ("unbound type parameter: '"^a))
      | TFun (_loc, _, _) -> <:expr<pp_fun>>
    and go_record fields =
      let go (_loc, fname, _, typ) =
        <:expr<
          $uid:if Hashtbl.mem omit fname then "Optional" else "Required"$
          ($qualify fname$,
           pp_precompose (fun x -> x.$lid:fname$) $go_type_expr typ$)>>
      in
      <:expr<fun ?(prec=0) fmt ->
          pp_record [| $fold_semi _loc (List.map go fields)$ |] fmt>>
    and go_variant variants =
      let f (_loc, ctor, args) =
        let ctor_str = <:expr<$uid:if Hashtbl.mem omit ctor
                                   then "Optional" else "Required"$
                              $qualify ctor$>>
        in
        match args with
        | [] -> <:match_case<$uid:ctor$ -> pp_string_noquote fmt $`str:ctor$>>
        | [t] -> <:match_case<$uid:ctor$ x ->
                    pp_ctor $ctor_str$
                            $go_type_expr t$
                            ~prec fmt x>>
        | ts ->
          let xs = mapi (fun i _ -> "x" ^ string_of_int i) args in
          let pat = fold_app_patt <:patt<$uid:ctor$>>
                      (List.map (fun x -> <:patt<$lid:x$>>) xs)
          and pps = fold_semi _loc
                       (List.map2
                         (fun x t -> <:expr<pp_fixarg $go_type_expr t$
                                            $lid:x$>>)
                         xs ts)
          in
          <:match_case<
            $pat$ ->
            pp_ctor $ctor_str$
              (fun ?(prec=0) fmt () ->
                pp_seq "(@[<hv>" "," "@])" fmt (Fstream.of_array [| $pps$ |]))
              ~prec fmt ()
          >>
      in
      <:expr<fun ?(prec=0) fmt -> function
        $fold_or_match_case _loc (List.map f variants)$>>
    in
    match rhs with
    | TRecord (_loc, fields) -> go_record fields
    | TVariant (_loc, variants) -> go_variant variants
    | TAlias rhs -> go_type_expr rhs

  (* Pre-defined aliases.  *)
  let pretty_predef_aliases = Hashtbl.create 10
  let _ =
    Hashtbl.add pretty_predef_aliases ["Lazy"; "t"] ["lazy_t"]

  (* Generate a family of pretty-printers for a set of types declared
     simultaneously.  The types may refer to each other, so the whole family is
     defined in a single let rec.  *)
  let pretty_of_type_decl_simpl tmp_aliases prefix omit (TDecl (_loc, types)) =
    let pretty_binding (_loc, params, tctor, rhs) =
      (* Load predefined aliases.  *)
      let aliases = Hashtbl.copy pretty_predef_aliases in
      Hashtbl.iter (Hashtbl.add aliases) tmp_aliases;
      let tctor =
        try
          let parts = Hashtbl.find aliases [tctor]
          in String.concat "_" parts
        with Not_found -> tctor in
      let pps = mapi (fun i x -> "pp"^string_of_int i) params in
      let wrap f =
        let complete_pp = fold_app <:expr<$lid:"pp_"^tctor$>>
                            (List.map (fun x -> <:expr<$lid:x$>>) pps) in
        fold_fun _loc pps <:expr<fun x -> $lid:f$ $complete_pp$ x>>
      in
      let pp_body = pretty_body _loc aliases prefix omit params rhs in
      let pp_impl =
        match pps with
        | [] -> (* eta-expand to avoid compiler error *)
                <:expr<fun ?(prec=0) fmt -> $pp_body$ ~prec fmt>>
        | pps -> fold_fun _loc pps pp_body
      in
      <:binding<
        $lid:"pp_"^tctor$ : $pretty_type _loc params tctor$ = $pp_impl$
        and $lid:"show_"^tctor$ = $wrap "show_of_pp"$
        and $lid:"dump_"^tctor$ = $wrap "dump_of_pp"$
        and $lid:"display_"^tctor$ = $wrap "display_of_pp"$
        and $lid:"ppout_"^tctor$ = $wrap "ppout_of_pp"$
        and $lid:"pperr_"^tctor$ = $wrap "pperr_of_pp"$
      >>
    in
    <:str_item<let rec
               $fold_and_binding _loc (List.map pretty_binding types)$>>

  (* Code generation options *)

  type pretty_options =
    { pretty_optional: (Loc.t * string list) option;
      pretty_prefix: (Loc.t * ident) option;
      pretty_aliases: (Loc.t * (string list * string list) list) option;
      pretty_erase_types: (Loc.t * bool) option; }

  let pretty_options_merge opt1 opt2 =
    let merge_field field_name field =
      match field opt1, field opt2 with
      | None, None -> None
      | Some _ as l, None -> l
      | None, (Some _ as r) -> r
      | Some (_loc1, _), Some (_loc2, _) ->
        fail_at (Loc.merge _loc1 _loc2) ("duplicate ~"^field_name^" options")
    in
    { pretty_optional = merge_field "optional" (fun x -> x.pretty_optional);
      pretty_prefix = merge_field "prefix" (fun x -> x.pretty_prefix);
      pretty_aliases = merge_field "alias" (fun x -> x.pretty_aliases);
      pretty_erase_types = merge_field "alias" (fun x -> x.pretty_erase_types);
    }

  let pretty_options_get f default opt =
    match opt with
    | None -> default
    | Some (_, x) -> f x

  let pretty_empty_options =
    { pretty_optional = None; pretty_prefix = None; pretty_aliases = None;
      pretty_erase_types = None; }

  EXTEND Gram
    GLOBAL: str_item module_longident type_longident;

    ident:
      [ [ `LIDENT x -> x
        | `UIDENT x -> x ]
      ];

    pretty_erase_typedefs:
      [ [ `LABEL "erase_typedefs"; "true" ->
          { pretty_empty_options with
            pretty_erase_types = Some (_loc, true) } ]
      ];

    pretty_aliases:
      [ [ `LABEL "rename"; "(";
          aliases = LIST1 [ x = type_longident; "to"; y = type_longident ->
                            (string_list_of_ident x, string_list_of_ident y) ]
                    SEP ",";
          ")" ->
          { pretty_empty_options with
            pretty_aliases = Some (_loc, aliases) } ]
      ];

    pretty_optional:
      [ [ `LABEL "optional"; "("; fs = LIST1 [ x = ident -> x ] SEP ","; ")" ->
          { pretty_empty_options with
            pretty_optional = Some (_loc, fs) } ]
      | [ `LABEL "optional"; x = ident ->
          { pretty_empty_options with
            pretty_optional = Some (_loc, [x]) } ]
      ];

    pretty_prefix:
      [ [ `LABEL "prefix"; x = module_longident ->
          { pretty_empty_options with
            pretty_prefix = Some (_loc, x) }
        ]
      ];

    deriving_spec:
      [ [ "("; `LIDENT "pretty";
          opts = LIST0 [ pretty_prefix | pretty_optional | pretty_aliases
                       | pretty_erase_typedefs ];
          ")" ->
          let opts = List.fold_left pretty_options_merge
                       pretty_empty_options opts
          in
          let optional =
            pretty_options_get (fun x -> x) [] opts.pretty_optional
          and prefix =
            pretty_options_get string_list_of_ident [] opts.pretty_prefix
          and aliases =
            pretty_options_get hashtbl_of_assocs (Hashtbl.create 0)
              opts.pretty_aliases
          and erase_types =
            pretty_options_get (fun x -> x) false opts.pretty_erase_types
          in
          fun td ->
            let pretty =
              pretty_of_type_decl_simpl aliases prefix
                (hashtbl_of_list optional) (type_decl_simpl_of_ctyp td)
            in
            if erase_types then pretty
            else <:str_item<type $td$
                            $pretty$>>
        | `LIDENT "pretty" ->
          fun td ->
            let pretty = pretty_of_type_decl_simpl
                           (Hashtbl.create 0)
                           []
                           (Hashtbl.create 0)
                           (type_decl_simpl_of_ctyp td)
            in <:str_item<type $td$ $pretty$>>
        ]
      ];

    str_item: LEVEL "top"
      [ "top"
        ["type"; td = type_declaration; "deriving"; spec = deriving_spec ->
         spec td
        ]
      ];

  END

end

let module M = Camlp4.Register.OCamlSyntaxExtension(Id)(Make) in ()
