(* A camlp4 extension for generating pretty printers and/or expr_of functions
for data types.  Requires OCaml 3.12 or newer.

[Quick Example]

The code

  open Pprint
  open Expr_of
  type 'a foo = Foo of ('a, int) bar
  deriving (pretty, expr_of)

defines the type foo and the functions

  pp_foo : 'a pp -> 'a foo pp                   (* pretty-printer for foo *)
  expr_of_foo : ('a -> expr) -> 'a foo -> expr  (* reifier for foo *)

where pp is the pretty-printer type defined in pprint.mli (which see; the other
pretty-printers defined there like show are also derived as show_foo).  The
type foo takes a type parameter, so the pp_foo and expr_of_foo take one
parameter each as well.  In general, a type with n parameters will have a
pretty-printer/reifier with n arguments.

If you want to derive the functions for a type that is already defined
somewhere else, then you can use the external_types keyword, like:

  open Pprint
  open Expr_of
  type 'a foo = Foo of ('a, int) bar
  deriving external_types (pretty, expr_of)

This defines the functions pp_foo and expr_of_foo as before, but the type
definition is erased from the preprocessed program source, so if foo is defined
in a different module, then the type is not re-defined.

"open Pprint" and "open Expr_of" are mandatory.  If the definition of foo
refers to another type baz, for example, the generated pretty-printer for foo
will call a pretty-printer named pp_baz.  The user can either write pp_baz by
hand or derive it with another deriving clause at the definition of type baz,
but in any case s/he must ensure that pp_baz exists and is visible.  Pprint
supplies these functions for base types like int, bool, etc., and it also
supplies some auxiliary functions that this generator relies on.  The same goes
for expr_of.

[More Detailed Semantics]

type t1 = ...
and  t2 = ...
deriving pretty

type t1 = ...
and  t2 = ...
deriving expr_of

type t1 = ...
and  t2 = ...
deriving (pretty, expr_of)

  These constructs derive just the pretty-printer, just expr_of, and both the
  pretty-printer and expr_of, respectively, for the types that are defined in
  the type definitions that appear just before the "deriving" keyword.

  If the type definition refers to a type t that is defined elsewhere, then
  the generated pretty-printer uses a function name pp_t to print that type,
  and the generated expr_of uses expr_of_t to reify that type.  It is up to the
  user to make sure that functions of those names are in scope.  If the type t
  has a module qualification, like Foo.t, then the inferred names are Foo.pp_t
  and Foo.expr_of_t, respectively.  As an exception, Lazy.t is handled by
  pp_lazy_t and expr_of_lazy_t (see the ~alias option below).

type t1 = ...
and  t2 = ...
deriving external_types pretty

type t1 = ...
and  t2 = ...
deriving external_types expr_of

type t1 = ...
and  t2 = ...
deriving external_types (pretty, expr_of)

  "deriving external_types" defines the pretty-printers and/or expr_of but
  eliminates the type definitions from the preprocessed source code.  This is
  useful for deriving from types that are defined in a different module, and
  you can't or don't want to modify that module.

  For example, option is defined in the Pervasives module, and you can't modify
  it.  But you can do

    (* In a module that you control:  *)
    type 'a option = None | Some of 'a
    deriving external_types (pretty, expr_of)

  and this will define just the pretty-printer and expr_of without re-defining
  the option type.  Note that pp_option and expr_of_option are already
  available in Pprint and Expr_of, respectively, so there's no point in doing
  this for option; you should do this on other types not covered by Pprint or
  Expr_of.

[Options]

The deriving commands pretty and expr_of can be given options, with syntax that
looks like optional function arguments.  Their syntax and semantics follow.

~optional:<constructor / record field / module name>
~optional:(<comma-separated list of constructor / record field / module name>)

  This is only available for pretty.  If more than one name is supplied, they
  must be separated by commas and the list of names must be parenthesized.
  Constructor names, records, or module names declared ~optional will be
  omitted from the output unless !Pprint.read_write_invariance is true.  For
  example,

      type foo = { important : int;
                   unimportant : int; }
      deriving (pretty ~optional:unimportant)

      let _ = pp_foo Format.std_formatter { important = 42; unimportant = 57 }

  prints

      {important = 42}

  if Pprint.read_write_invariance is false, but the same code prints

      {important = 42; unimportant = 57}

  if Pprint.read_write_invariance is true.  This is useful for hiding internal
  data that is not interesting to the user and keeping the output concise. For
  optional constructors, the constructor name is omitted but its arguments are
  not; for example,

      type positive_int = Positive of int
      deriving (pretty ~optional:Positive)

      let _ = pp_foo Format.std_formatter (Positive 5)

  prints "5" if Pprint.read_write_invariance is false, and the same code prints
  "Positive 5" if Pprint.read_write_invariance is true.

~prefix:<modname>

  The argument must not be parenthesized.  <modname> is prepended to every
  constructor or field name, so that the pretty-printed result or the output of
  expr_of is useful outside of the current module.  Don't confuse this with
  ~opening below, which prepends a module name to values and names used in the
  generated functions; by contrast, ~prefix attaches a module name to the names
  appearing the outputs of the generated functions.

  For example, the code

      (* In foo.ml *)
      type x = X
      deriving (pretty, expr_of)

      (* In bar.ml *)
      type x = X

      (* In baz.ml *)
      let s = Foo.show_x Foo.X
      let e = Foo.expr_of_x Foo.X

  binds s to "X" and e to <:expr<X>>.  But does this "X" mean Foo.X or Bar.X?
  You probably mean Foo.X, but this is not reflected in the values.  To avoid
  this ambiguity, we can write

      (* In foo.ml *)
      type x = X
      deriving (pretty ~prefix:Foo, expr_of ~prefix:Foo)

      (* In bar.ml *)
      type x = X

      (* In baz.ml *)
      let s = Foo.show_x Foo.X
      let e = Foo.expr_of_x Foo.X

  which binds s to "Foo.X" and e to <:expr<Foo.X>>.

  In pretty, the prefix is subject to hiding by ~optional, so for example

      type x = X
      deriving (pretty ~prefix:Foo ~optional:Foo)

      let _ = pp_x Format.std_formatter X

  prints "X" if Pprint.read_write_invariance is false but prints "Foo.X" if
  Pprint.read_write_invariance is true.

~opening:<modname>

  The argument must not be parenthesized.  The derived functions are defined
  with "let open <modname>" around them.  Don't confuse this with ~prefix,
  which attaches a module name to names in the outputs of the generated
  functions, whereas ~opening attaches a module name to the values and types
  used in the generated functions.

  This is useful in combination with external_types.  For example,

    (* In foo.ml *)
    type foo = Foo

    (* In bar.ml *)
    type foo = Foo
    deriving external_types (pretty ~opening:Foo)

  Without the ~opening:Foo, the code in bar.ml would produce an error because
  it would generate code that looks like

    (* In bar.ml *)
    let rec pp_foo : foo pp (* error: foo is not in scope *)
     = ...

  With ~opening:Foo, this is changed to

    include struct
      open Foo
      let rec pp_foo : foo pp (* OK *)
       = ...
    end

  (BTW, the type annotation is needed in the first place to ensure the right
   degree of polymorphism.)

~alias:(<name> = <name>, <name> = <name>, ...)

  This option gives a different name from which to guess the name of the
  pretty-printer or expr_of function of a type.  For example,

      type 'a foo = Foo of 'a ThirdPartyModule.t
      deriving pretty

  generates a pretty-printer that refers to ThirdPartyModule.pp_t which may not
  be supplied by whoever wrote ThirdPartyModule.  In that case, you can use
  ~alias to make deriving pretend that ThirdPartyModule.t has a different name:

      let pp_tpm_t pp_'a = (* pretty-printer for 'a ThirdPartyModule.t *)
      type foo = Foo of 'a ThirdPartyModule.t
      deriving (pretty ~alias:(ThirdPartyModule.t = tpm_t))

  Then pp_foo will use pp_tpm_t to print 'a ThirdPartyModule.t instead of the
  non-existent ThirdPartyModule.pp_t function.

  The alias Lazy.t = lazy_t is loaded automatically for every deriving clause.

[Syntax]

The complete syntax for deriving is as follows, using pseudo-EBNF with the
usual convention that [tok] denotes an optional token, * is the Kleene star
(i.e. repetition for 0 or more times), and | denotes alternatives.  Parentheses
do not have special meanings; they denote parentheses.  LIST <sep> <foo> means
a <sep>-separated list of 1 or more <foo>'s.

<deriving_clause> ::=
  <type def>
  deriving [external_types] <deriving_commands>
  [;;]

<deriving_commands> ::=
  pretty | expr_of | ( LIST , <deriving_command_with_opt> )

<deriving_command_with_opt> ::=
  pretty <pretty_option>* | expr_of <expr_of_option>*

<pretty_option> ::=
  ~alias:( LIST , <alias> ) | ~prefix:<modname> | ~optional:<names>

<expr_of_option> ::=
  ~alias:( LIST , <alias> ) | ~prefix:<modname>

<alias> ::= <name> = <name>

<names> ::= <name> | ( LIST , <name> )

<modname> is any OCaml module name.  For instance Foo is OK (and it doesn't
matter if a module by that name exists), Foo.Bar is OK, but Foo.bar is not.

<name> is any valid OCaml identifier, capitalized or not.  It can contain
module qualifications like Foo.bar if needed.

<type def> is any valid type definition of 1 or more non-opaque types.  So
"type foo = int and 'a bar = 'a" is OK, but "type foo" (without a right-hand
side) is not OK.

*)

module Id : Camlp4.Sig.Id = struct
  let name = "pa_deriving"
  let version = "1.0"
end

module Make (Syntax : Camlp4.Sig.Camlp4Syntax) =
struct
  open Camlp4.Sig
  include Syntax
  open Ast

  module Camlp4aux = Camlp4aux.Make (Syntax)
  open Camlp4aux

  module M = Meta.Make (Meta.MetaLoc)
  module MExpr = M.Expr
  module MPatt = M.Patt

  (*
   * Auxiliary defs
   *)

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
    | <:ctyp<$_$ = $t$>> -> type_def_simpl_of_ctyp t
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
        let x = id_of_string_list _loc (map_last (fun y -> "pp_"^y) x) in
        fold_app <:expr<$id:x$>> (List.map go_type_expr args)
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
            pp_seq "(@[<hv>" "," "@])" fmt
               (Fstream.of_array [| $fold_sem _loc (mapi mk_pp xs)$ |])
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
          pp_record [| $fold_sem _loc (List.map go fields)$ |] fmt>>
    and go_variant variants =
      let f (_loc, ctor, args) =
        let ctor_str = <:expr<$uid:if Hashtbl.mem omit ctor
                                   then "Optional" else "Required"$
                              $qualify ctor$>>
        in
        match args with
        | [] -> <:match_case<$uid:ctor$ ->
                             pp_string_noquote fmt $qualify ctor$>>
        | [t] -> <:match_case<$uid:ctor$ x ->
                    pp_ctor $ctor_str$
                            $go_type_expr t$
                            ~prec fmt x>>
        | ts ->
          let xs = mapi (fun i _ -> "x" ^ string_of_int i) args in
          let pat = fold_app_patt <:patt<$uid:ctor$>>
                      (List.map (fun x -> <:patt<$lid:x$>>) xs)
          and pps = fold_sem _loc
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
  let predef_aliases = Hashtbl.create 10
  let _ =
    Hashtbl.add predef_aliases ["Lazy"; "t"] ["lazy_t"]

  (* Generate a family of pretty-printers for a set of types declared
     simultaneously.  The types may refer to each other, so the whole family is
     defined in a single let rec.  *)
  let pretty_of_type_decl_simpl tmp_aliases prefix omit (TDecl (_loc, types)) =
    let pretty_binding (_loc, params, tctor, rhs) =
      (* Load predefined aliases.  *)
      let aliases = Hashtbl.copy predef_aliases in
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
               $biAnd_of_list (List.map pretty_binding types)$>>

  (* Reification (expr_of) derivation  *)

  (* Generate a type for the reifier of the given type.  This is
     necessary since the reifier for e.g.
       type 'a foo = ...
       and bar = int foo
     requires polymorphic recursion.  This is the part that needs
     OCaml >= 3.12.  *)
  let expr_type _loc param_names tctor =
    let params = List.map (fun x -> <:ctyp<'$lid:x$>>) param_names in
    forall _loc param_names
      (fold_arr (List.map (fun x -> <:ctyp<$x$ -> Ast.expr>>) params)
         (TyArr (_loc, fold_app_ctyp <:ctyp<$lid:tctor$>> params,
                 <:ctyp<Ast.expr>>)))

  (* Generate a reifier's body.  The resulting expression has type 'a ->
     expr.  *)
  let expr_body _loc aliases prefix params rhs =
    (* Table from type parameter name to expr0, expr1, etc. *)
    let ex_table = Hashtbl.create 10 in
    List.iter (fun x -> Hashtbl.add ex_table x
                  <:expr<$lid:"expr"^string_of_int
                           (Hashtbl.length ex_table)$>>)
      params;
    let ghost = <:expr<Loc.ghost>> in
    let prefix = List.map (fun x -> <:ident<$uid:x$>>) prefix in
    let lqualify str = idAcc_of_list (prefix @ [ <:ident<$lid:str$>> ]) in
    let uqualify str = idAcc_of_list (prefix @ [ <:ident<$uid:str$>> ]) in
    let rec go_type_expr = function
      | TApp (_loc, <:ident<$id:x$>>, args) ->
        let x = string_list_of_ident x in
        let x = try Hashtbl.find aliases x with Not_found -> x in
        let x = id_of_string_list _loc (map_last (fun y -> "expr_of_"^y) x) in
        fold_app <:expr<$id:x$>> (List.map go_type_expr args)
      | TTup (_loc, []) -> fail_at _loc
                           "internal error in deriving: TTup (_,[])"
      | TTup (_loc, [x]) -> fail_at _loc
                           "internal error in deriving: TTup (_,[_])"
      | TTup (_loc, [x;y]) -> <:expr<expr_of_pair $go_type_expr x$
                                                  $go_type_expr y$>>
      | TTup (_loc, [x;y;z]) -> <:expr<expr_of_triple $go_type_expr x$
                                                      $go_type_expr y$
                                                      $go_type_expr z$>>
      | TTup (_loc, xs) ->
        let mk_pat i _ = <:patt<$lid:"x"^string_of_int i$>>
        and mk_expr i x = <:expr<$go_type_expr x$ $lid:"x"^string_of_int i$>>
        in
        <:expr<function
          | $fold_tuple_patt (mapi mk_pat xs)$ ->
            $meta_tuple _loc ghost (mapi mk_expr xs)$
        >>
      | TVar (_loc, a) ->
        (try Hashtbl.find ex_table a
         with Not_found -> fail_at _loc ("unbound type parameter: '"^a))
      | TFun (_loc, _, _) -> <:expr<pp_fun>>
    and go_record _loc fields =
      let go (_loc, fname, _, typ) =
        let qfname = MExpr.meta_ident _loc (lqualify fname) in
        <:expr< Ast.RbEq ($ghost$, $qfname$,
                          $go_type_expr typ$ x.$lid:fname$) >>
      in
      let binding = meta_fold_sem_rec_binding ghost (List.map go fields) in
      <:expr<fun x -> Ast.ExRec ($ghost$, $binding$, Ast.ExNil $ghost$) >>
    and go_variant variants =
      let f (_loc, ctor, args) =
        let qctor = MExpr.meta_expr _loc (meta_ident (uqualify ctor)) in
        match args with
        | [] -> <:match_case<$uid:ctor$ -> $qctor$>>
        | [t] -> <:match_case<$uid:ctor$ x ->
                              $(meta_app ghost qctor
                                (app (go_type_expr t) (meta_lid _loc "x")))$>>
        | ts ->
          let xs = mapi (fun i _ -> "x" ^ string_of_int i) args in
          let pat = fold_app_patt <:patt<$uid:ctor$>>
                      (List.map (fun x -> <:patt<$lid:x$>>) xs)
          and args = List.map2
                       (fun x t -> <:expr<$go_type_expr t$ $lid:x$>>)
                       xs ts
          in
          let body = meta_fold_app ghost (qctor::args) in
          <:match_case<$pat$ -> $body$>>
      in
      <:expr<function $fold_or_match_case _loc (List.map f variants)$>>
    in
    match rhs with
    | TRecord (_loc, fields) -> go_record _loc fields
    | TVariant (_loc, variants) -> go_variant variants
    | TAlias rhs -> go_type_expr rhs

  let expr_of_of_type_decl_simpl _loc prefix tmp_aliases (TDecl (_loc, types))
    =
    let expr_binding (_loc, params, tctor, rhs) =
      (* Load predefined aliases.  *)
      let aliases = Hashtbl.copy predef_aliases in
      Hashtbl.iter (Hashtbl.add aliases) tmp_aliases;
      let tctor =
        try
          let parts = Hashtbl.find aliases [tctor]
          in String.concat "_" parts
        with Not_found -> tctor in
      let exprs = mapi (fun i x -> "expr"^string_of_int i) params in
      let expr_body = expr_body _loc aliases prefix params rhs in
      let expr_impl =
        match exprs with
        | [] -> (* eta-expand to avoid compiler error *)
                <:expr<fun x -> $expr_body$ x>>
        | exprs -> fold_fun _loc exprs expr_body
      in
      <:binding<
        $lid:"expr_of_"^tctor$ : $expr_type _loc params tctor$ = $expr_impl$
      >>
    in
    <:str_item<let rec
               $fold_and_binding (List.map expr_binding types)$>>

  let wrap_opening opening str_item =
    match opening with
    | None -> str_item
    | Some (_loc, modname) ->
      <:str_item<include struct
                   open $id:modname$
                   $str_item$
                 end>>

  (* Code generation options *)

  type options =
    { opt_optional: (Loc.t * string list) option;
      opt_prefix: (Loc.t * ident) option;
      opt_aliases: (Loc.t * (string list * string list) list) option;
      opt_opening: (Loc.t * ident) option;
    }

  let options_merge opt1 opt2 =
    let merge_field field_name field =
      match field opt1, field opt2 with
      | None, None -> None
      | Some _ as l, None -> l
      | None, (Some _ as r) -> r
      | Some (_loc1, _), Some (_loc2, _) ->
        fail_at (Loc.merge _loc1 _loc2) ("duplicate ~"^field_name^" options")
    in
    { opt_optional = merge_field "optional" (fun x -> x.opt_optional);
      opt_prefix = merge_field "prefix" (fun x -> x.opt_prefix);
      opt_aliases = merge_field "alias" (fun x -> x.opt_aliases);
      opt_opening = merge_field "opening" (fun x -> x.opt_opening);
    }

  let options_get f default opt =
    match opt with
    | None -> default
    | Some (_, x) -> f x

  let options_empty =
    { opt_optional = None; opt_prefix = None; opt_aliases = None;
      opt_opening = None;
    }

  (* Entry points for derivations  *)

  let derive_pretty _loc opts =
    let optional = options_get (fun x -> x) [] opts.opt_optional
    and prefix = options_get string_list_of_ident [] opts.opt_prefix
    and aliases = options_get hashtbl_of_assocs (Hashtbl.create 0)
        opts.opt_aliases
    in
    fun td ->
      let pretty =
        pretty_of_type_decl_simpl aliases prefix
          (hashtbl_of_list optional) (type_decl_simpl_of_ctyp td)
      in
      wrap_opening opts.opt_opening
        <:str_item<$pretty$>>

  let derive_expr_of _loc opts =
    let prefix = options_get string_list_of_ident [] opts.opt_prefix
    and aliases = options_get hashtbl_of_assocs (Hashtbl.create 0)
        opts.opt_aliases
    in
    fun td ->
      wrap_opening opts.opt_opening
        (expr_of_of_type_decl_simpl _loc prefix aliases
           (type_decl_simpl_of_ctyp td))

  EXTEND Gram
    GLOBAL: str_item module_longident type_longident;

    ident:
      [ [ `LIDENT x -> x
        | `UIDENT x -> x ]
      ];

    label_aliases:
      [ [ `LABEL "alias" -> () | `LABEL "aliases" -> () ] ];

    opt_aliases:
      [ [ label_aliases; "(";
          aliases = LIST1 [ x = type_longident; "="; y = type_longident ->
                            (string_list_of_ident x, string_list_of_ident y) ]
                    SEP ",";
          ")" ->
          { options_empty with
            opt_aliases = Some (_loc, aliases) } ]
      ];

    opt_optional:
      [ [ `LABEL "optional"; "("; fs = LIST1 [ x = ident -> x ] SEP ","; ")" ->
          { options_empty with
            opt_optional = Some (_loc, fs) } ]
      | [ `LABEL "optional"; x = ident ->
          { options_empty with
            opt_optional = Some (_loc, [x]) } ]
      ];

    opt_prefix:
      [ [ `LABEL "prefix"; x = module_longident ->
          { options_empty with
            opt_prefix = Some (_loc, x) }
        ]
      ];

    opt_opening:
      [ [ `LABEL "opening"; x = module_longident ->
          { options_empty with
            opt_opening = Some (_loc, x) }
        ]
      ];

    bare_pretty_spec:
      [ [ `LIDENT "pretty" -> derive_pretty _loc options_empty ]
      ];

    pretty_spec:
      [ [ `LIDENT "pretty";
          opts = LIST0 [ opt_prefix | opt_optional | opt_aliases
                       | opt_opening ]
          ->
          derive_pretty _loc (List.fold_left options_merge options_empty opts)
        ]
      ];

    bare_expr_of_spec:
      [ [ `LIDENT "expr_of" -> derive_expr_of _loc options_empty ]
      ];

    expr_of_spec:
      [ [ `LIDENT "expr_of"; opts = LIST0 [ opt_prefix | opt_aliases
                                          | opt_opening ]
          ->
          derive_expr_of _loc (List.fold_left options_merge options_empty opts)
        ]
      ];

    deriving_spec:
      [ [ bare_pretty_spec | bare_expr_of_spec
        ]
      | [ "("; specs = LIST1 [ pretty_spec | expr_of_spec ] SEP ","; ")" ->
          fun td -> fold_str_item _loc (List.map (fun f -> f td) specs)
        ]
      ];

    str_item: LEVEL "top"
      [ "top"
        ["type"; td = type_declaration; "deriving";
         erase = OPT [ `LIDENT "external_types" ]; spec = deriving_spec ->
         if erase = None then <:str_item<type $td$ $spec td$>>
         else spec td
        ]
      ];

  END

end

let module M = Camlp4.Register.OCamlSyntaxExtension(Id)(Make) in ()
