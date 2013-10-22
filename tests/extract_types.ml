(* camlp4orf *)
(* Outside quotes -> original syntax
   Inside quotes -> revised syntax *)

open Camlp4
open PreCast
open Syntax
open Ast
module RParser = Camlp4OCamlRevisedParser.Make(Syntax)
module Parser = Camlp4OCamlParser.Make(Syntax)
module Camlp4aux = Camlp4aux.Make(Syntax)

let ocaml_printer =
  Camlp4aux.make_printer ~comments:false ~semisep:"" ()

let print_str_item = ocaml_printer#print_str_item
let string_of_str_item = ocaml_printer#string_of_str_item
let print_sig_item = ocaml_printer#print_sig_item
let string_of_sig_item = ocaml_printer#string_of_sig_item
let print_ctyp = ocaml_printer#print_ctyp
let string_of_ctyp = ocaml_printer#string_of_ctyp

let type_templ = ref ""
and exc_templ = ref ""
and exceptions = ref false
and file = ref ""
and source = ref ""
and modname = ref ""
and debug = ref false
and header = ref ""
and footer = ref ""
and excluded = ref []
and debug_exclude = ref false

let die msg =
  Printf.fprintf stderr "Error: %s\n" msg;
  exit 1

let warn msg =
  Printf.fprintf stderr "Warning: %s\n" msg;
  flush stderr

let dbg msg = if !debug then (Printf.fprintf stderr "Debug: %s\n" msg;
                              flush stderr)

let set_string_check_dup opt dest s =
  if !dest = "" then dest := s
  else die (Printf.sprintf "%s specified more than once (%s and %s)"
           opt !dest s)
let string_opt opt dest doc =
  (opt, Arg.String (set_string_check_dup opt dest), doc)

let process_argv () =
  Arg.parse
    [("--exceptions", Arg.Set exceptions,
      "extract exceptions too");
     string_opt "--replace" type_templ
       "replace each type declaration by a formatted string";
     ("--replace-exc", Arg.String (fun s ->
          set_string_check_dup "--replace-exc" exc_templ s;
          exceptions := true),
      "replace each exception declaration by a formatted string; \
       implies --exceptions");
     string_opt "--source" source
       "set the name of the source file, if different from input file name";
     ("--debug", Arg.Set debug, "print debugging information to stderr");
     ("--header", Arg.Set_string header,
      "string to print before all other outputs");
     ("--footer", Arg.Set_string footer,
      "string to print after all other outputs");
     ("--exclude",
      Arg.String (fun s -> excluded := Str.regexp s :: !excluded),
      "exclude types/exceptions matching a regexp");
     ("--debug-exclude",
      Arg.Set debug_exclude,
      "print to stderr which types are excluded and which are kept");
    ]
    (set_string_check_dup "input file" file)
    "\
A utility that extracts only the type and exception definitions from a given
.ml(i) file and prints them to stdout.  Used to extract pretty-printers and
expr_of functions from a module without modifying it.

The flags --replace and --replace-exc take a string as argument, and each type
or exception declaration is replaced by that string.  The string can contain
the following format directives:

    %M The full module name in which the type/exception appears.
    %m The (sub)module name, not including the name of the .ml(i) file.
    %T The type or exception declaration, with a newline before and after.
    %t The type or exception declaration with no newlines.
    %% The character '%'.

The following escape characters are also recognized:

    \\n
    \\t
    \\\\

Both --replace and --replace-exc default to %T if omitted.

This utility doesn't do any camlp4 translation on the .ml(i) file it reads, so
the file must not use syntax extensions; if it does use extensions, the
file should be pre-processed with camlp4 and printed (in human-readable
format) before feeding it to this utility.  The --source option may come in
handy in that case.

Unwanted types can be filtered out by --exclude <regexp>.  The name of each
type or exception, qualified with submodule names if it's declared inside a
submodule, is matched against <regexp> and omitted from the output if it
matches.  When the option is given more than once, a type or exception is
excluded if it matches at least one of the <regexp>s.
";
  if !file = "" then die "no input file";
  if !source = "" then source := !file;
  if !exc_templ = "" then exc_templ := "%T";
  if !type_templ = "" then type_templ := "%T";
  let modname_regexp = Str.regexp "^[A-Z][a-zA-Z0-9_]*$" in
  let filename_regexp = Str.regexp "^\\([a-zA-Z][a-zA-Z0-9_]*\\)\\.mli?$" in
  (* Extract module name from file name.  *)
  if Str.string_match filename_regexp (Filename.basename !file) 0 then
    begin
      modname := Str.matched_group 1 (Filename.basename !file);
      !modname.[0] <- Char.uppercase !modname.[0];
      dbg (Printf.sprintf "module name = (%s), file name = (%s)"
             !modname !file);
      if not (Str.string_match modname_regexp !modname 0)
      then die (Printf.sprintf "internal error: file name %s appears legit \
                                but module name %s is broken" !file !modname)
    end
  else die (Printf.sprintf "File name %s doesn't give a valid module name."
              !file)

let exclude submod tctor =
  let type_name = String.concat "." (submod @ [tctor]) in
  let match_against r =
    try let _ = Str.search_forward r type_name 0 in
        true
    with Not_found -> false
  in
  let ret = List.exists match_against !excluded in
  if !debug_exclude then
    Printf.fprintf stderr "%s type %s\n"
      (if ret then "excluded" else "kept")
      type_name;
  ret

let rec trim_tydcl submod td =
  match td with
  | TyDcl (_, _, _, (TyNil _ as nil), _) ->
    (* <<type foo>> with no rhs; ignore *)
    nil
  | TyDcl (_loc, tctor, _, _, _) when exclude submod tctor -> TyNil _loc
  | TyDcl (_, _, _, _, _) -> td
  | TyAnd (_loc, l, r) ->
    (match trim_tydcl submod l, trim_tydcl submod r with
     | TyNil _, TyNil _ -> TyNil _loc
     | TyNil _, td -> td
     | td, TyNil _ -> td
     | l, r -> TyAnd (_loc, l, r))
  | _ -> die (Printf.sprintf "unrecognized type definition %s"
                (string_of_ctyp td))

let trim_exc submod exc =
  match exc with
  | TyOf (_, TyId (_, IdUid (_, ctor)), _) | TyId (_, IdUid (_, ctor)) ->
    if !exceptions && not (exclude submod ctor)
    then exc
    else TyNil (loc_of_ctyp exc)
  | _ -> die (Printf.sprintf "unrecognized exception definition %s"
                (string_of_ctyp exc))

(* The main type extraction engine.  Extracts a list of type and exception
   declarations along with the list of modules encountered along the way.  *)
let sig_extract_types =
  object (self)
    inherit Ast.fold as super

    val submod = []                     (* current submodule *)
    val acc = []

    method acc = List.rev acc

    method sig_item s =
      let _loc = loc_of_sig_item s in
      match s with
      | <:sig_item<type $td$>> ->
        (match trim_tydcl submod td with
         | TyNil _ -> self
         | td' -> {< acc = (submod, <:sig_item<type $td'$>>)::acc >})
      | <:sig_item<exception $exc$>> as s ->
        (match trim_exc submod exc with
         | TyNil _ -> self
         | exc' -> {< acc = (submod, s)::acc >})
      | <:sig_item<module $modname$ : sig $signature$ end>> ->
        let from_submod =
          ({< submod = submod@[modname] >}#sig_item signature)#acc
        in {< acc = List.rev_append from_submod acc >}
      | <:sig_item<class $_$>> ->
        warn (Printf.sprintf "%s: classes are not supported, ignoring..."
               (Loc.to_string _loc));
        self
      | s -> super#sig_item s
  end

let str_extract_types =
  object (self)
    inherit Ast.fold as super

    val submod = []                     (* current submodule *)
    val acc = []

    method acc = List.rev acc

    method str_item s =
      let _loc = loc_of_str_item s in
      match s with
      | <:str_item<type $td$>> ->
        (match trim_tydcl (List.rev submod) td with
         | TyNil _ -> self
         | td' -> {< acc = (List.rev submod, <:str_item<type $td'$>>)::acc >})
      | <:str_item<exception $exc$>> as s ->
        (match trim_exc submod exc with
         | TyNil _ -> self
         | exc' -> {< acc = (submod, s)::acc >})
      | <:str_item<module $modname$ = struct $structure$ end>> ->
        let from_submod =
          ({< submod = submod@[modname] >}#str_item structure)#acc
        in {< acc = List.rev_append from_submod acc >}
      | <:str_item<module $modname$ : $_$ = struct $structure$ end>> ->
        let from_submod =
          ({< submod = submod@[modname] >}#str_item structure)#acc
        in {< acc = List.rev_append from_submod acc >}
      | <:str_item<class $_$>> ->
        warn (Printf.sprintf "%s: classes are not supported, ignoring..."
               (Loc.to_string _loc));
        self
      | s -> super#str_item s
  end

let with_infile file thunk =
  match file with
  | "-" -> thunk stdin
  | _ ->
    let infile = try open_in file with Sys_error msg -> die msg in
    try
      let ret = thunk infile in
      close_in infile;
      ret
    with exn -> close_in infile; raise exn

let print_header_footer str =
  let rec parse =
    parser
    | [< ''\\'; str >] -> escape str
    | [< 'c; str >] -> print_char c; parse str
    | [< >] -> print_char '\n'
  and escape =
    parser
    | [< ''n'; str >] -> print_char '\n'; parse str
    | [< ''t'; str >] -> print_char '\t'; parse str
    | [< ''\\'; str >] -> print_char '\\'; parse str
    | [< 'c; _ >] -> die (Printf.sprintf
                            "unrecognized escape character \\%c in %s"
                            c str)
  in
  if str <> "" then parse (Stream.of_string str)

let print_tdecl print_item full_modname submodname item templ =
  let rec parse =
    parser
    | [< ''%'; str >] -> directive str
    | [< ''\\'; str >] -> escape str
    | [< 'c; str >] -> print_char c; parse str
    | [< >] -> ()
  and escape =
    parser
    | [< ''n'; str >] -> print_char '\n'; parse str
    | [< ''t'; str >] -> print_char '\t'; parse str
    | [< ''\\'; str >] -> print_char '\\'; parse str
    | [< 'c; _ >] -> die (Printf.sprintf
                            "unrecognized escape character \\%c in %s"
                            c templ)
  and directive =
    parser
    | [< ''T'; str >] -> print_char '\n';
                         print_item item;
                         print_char '\n';
                         parse str
    | [< ''t'; str >] -> print_item item; parse str
    | [< ''M'; str >] -> print_string full_modname; parse str
    | [< ''m'; str >] -> print_string submodname; parse str
    | [< ''%'; str >] -> print_char '%'; parse str
    | [< 'c; _ >] -> die (Printf.sprintf
                            "unrecognized %%-directive %%%c in %s" c templ)
    | [< >] -> die ("string ends with a %: " ^ templ)
  in parse (Stream.of_string templ)

let print_sig (submodname, item) =
  let full_modname = String.concat "." (!modname::submodname)
  and submodname = String.concat "." submodname in
  let print_tdecl = print_tdecl print_sig_item full_modname submodname in
  match item with
  | <:sig_item<type $_$>> -> print_tdecl item !type_templ
  | <:sig_item<exception $_$>> -> print_tdecl item !exc_templ
  | _ -> die ("internal error: extracted something other than \
               type or exception: " ^ string_of_sig_item item)

let print_str (submodname, item) =
  let full_modname = String.concat "." (!modname::submodname)
  and submodname = String.concat "." submodname in
  let print_tdecl = print_tdecl print_str_item full_modname submodname in
  match item with
  | <:str_item<type $_$>> -> print_tdecl item !type_templ
  | <:str_item<exception $_$>> -> print_tdecl item !exc_templ
  | _ -> die ("internal error: extracted something other than \
               type or exception: " ^ string_of_str_item item)

let parse parser_fcn srcname =
  try
    let loc = Loc.of_tuple (srcname, 1, 0, 0, 1, 0, 0, false) in
    with_infile srcname (fun infile ->
        parser_fcn loc (Stream.of_channel infile))
  with Loc.Exc_located (loc, Stream.Error msg) ->
    Printf.fprintf stderr "%s:\nError: %s\n" (Loc.to_string loc) msg;
    failwith "Input file could not be parsed."

let main () =
  (* print_interf seems to always fail with "No interface printer".  I've tried
     module Printer = Camlp4OCamlPrinter, I've also tried calling
     Camlp4.Register.enable_ocaml_printer () directly, and neither did anything
     at all.  If we go directly to the Printers.OCaml module it seems to
     work.  *)
  process_argv ();
  if Filename.check_suffix !file ".mli" then
    let src = parse Parser.parse_interf !file in
    let skeleton = (sig_extract_types#sig_item src)#acc in
    print_header_footer !header;
    List.iter print_sig skeleton;
    print_header_footer !footer;
    dbg (Printf.sprintf "collected %d types." (List.length skeleton))
  else if Filename.check_suffix !file ".ml" then
    let src = parse Parser.parse_implem !file in
    let skeleton = (str_extract_types#str_item src)#acc in
    print_header_footer !header;
    List.iter print_str skeleton;
    print_header_footer !footer;
    dbg (Printf.sprintf "collected %d types." (List.length skeleton))
  else
    die (Printf.sprintf "unrecognized extension: %s" !file)
;;

try
  main ()
with Loc.Exc_located (loc, exn) ->
  raise exn
