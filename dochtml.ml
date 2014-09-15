(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

(**
 Custom tags for the ocamldoc comments:
    @cvode          link to Sundials CVODE documentation
    @cvodes         link to Sundials CVODES documentation
    @ida            link to Sundials IDA documentation
    @idas           link to Sundials IDAS documentation
    @kinsol         link to Sundials KINSOL documentation
 *)

let cvode_doc_root =
  ref "https://computation.llnl.gov/casc/sundials/documentation/cv_guide/"

let cvodes_doc_root =
  ref "https://computation.llnl.gov/casc/sundials/documentation/cvs_guide/"

let ida_doc_root =
  ref "https://computation.llnl.gov/casc/sundials/documentation/ida_guide/"

let idas_doc_root =
  ref "https://computation.llnl.gov/casc/sundials/documentation/idas_guide/"

let kinsol_doc_root =
  ref "https://computation.llnl.gov/casc/sundials/documentation/kin_guide/"

#if OCAML_3X == 1
class dochtml =
  object(self)
    inherit Odoc_html.html as super
#else
module Generator (G : Odoc_html.Html_generator) =
struct
  class html =
  object(self)
    inherit G.html as super
#endif

    val rex = Str.regexp "<\\([^#>]*\\)\\(#[^)]*\\)?> \\(.*\\)"

    method private split_text (t:Odoc_info.text) =
      let s = Odoc_info.text_string_of_text t in
      if not (Str.string_match rex s 0) then
        failwith "Bad parse!"
      else
        let page = Str.matched_group 1 s
        and anchor = try
              Str.matched_group 2 s
            with Not_found -> ""
        and title = Str.matched_group 3 s
        in
      (page, anchor, title)

    method private html_of_cvode t =
      let (page, anchor, title) = self#split_text t in
      Printf.sprintf
        "<div class=\"cvode\"><small>See sundials: <a href=\"%s%s.html%s\">%s</a></small></div>"
        !cvode_doc_root page anchor title

    method private html_of_cvodes t =
      let (page, anchor, title) = self#split_text t in
      Printf.sprintf
        "<div class=\"cvodes\"><small>See sundials: <a href=\"%s%s.html%s\">%s</a></small></div>"
        !cvodes_doc_root page anchor title

    method private html_of_ida t =
      let (page, anchor, title) = self#split_text t in
      Printf.sprintf
        "<div class=\"ida\"><small>See sundials: <a href=\"%s%s.html%s\">%s</a></small></div>"
        !ida_doc_root page anchor title

    method private html_of_idas t =
      let (page, anchor, title) = self#split_text t in
      Printf.sprintf
        "<div class=\"idas\"><small>See sundials: <a href=\"%s%s.html%s\">%s</a></small></div>"
        !idas_doc_root page anchor title

    method private html_of_kinsol t =
      let (page, anchor, title) = self#split_text t in
      Printf.sprintf
        "<div class=\"kinsol\"><small>See sundials: <a href=\"%s%s.html%s\">%s</a></small></div>"
        !kinsol_doc_root page anchor title

    val divrex = Str.regexp " *\\(open\\|close\\) *\\(.*\\)"

    method private html_of_div text =
      let baddiv s = (Odoc_info.warning (Printf.sprintf
            "div must be followed by 'open' or 'close', not '%s'!" s); "") in
      match text with
      | [Odoc_info.Raw s] ->
          if not (Str.string_match divrex s 0) then baddiv s
          else begin
            match Str.matched_group 1 s with
            | "open" ->
                let attrs = try Str.matched_group 2 s with Not_found -> "" in
                Printf.sprintf "<div %s>" attrs
            | "close" -> "</div>"
            | s -> baddiv s
          end
      | _ -> failwith "div must be followed by plain text!"

    val mutable custom_functions =
      ([] : (string * (Odoc_info.text_element list -> string)) list)

    method private html_of_custom_text b tag text =
      try
        let f = List.assoc tag custom_functions in
        Buffer.add_string b (f text)
      with
        Not_found ->
          Odoc_info.warning (Odoc_messages.tag_not_handled tag)

    initializer
      tag_functions <- ("cvode",  self#html_of_cvode) :: tag_functions;
      tag_functions <- ("cvodes", self#html_of_cvodes) :: tag_functions;
      tag_functions <- ("ida",    self#html_of_ida) :: tag_functions;
      tag_functions <- ("idas",   self#html_of_idas) :: tag_functions;
      tag_functions <- ("kinsol", self#html_of_kinsol) :: tag_functions;

      custom_functions <- ("div", self#html_of_div) :: custom_functions

  end

#if OCAML_3X == 0
end
#endif

let option_cvode_doc_root =
  ("-cvode-doc-root", Arg.String (fun d -> cvode_doc_root := d), 
   "<dir>  specify the root url for the Sundials CVODE documentation.")
let option_cvodes_doc_root =
  ("-cvodes-doc-root", Arg.String (fun d -> cvodes_doc_root := d), 
   "<dir>  specify the root url for the Sundials CVODES documentation.")
let option_ida_doc_root =
  ("-ida-doc-root", Arg.String (fun d -> ida_doc_root := d), 
   "<dir>  specify the root url for the Sundials IDA documentation.")
let option_idas_doc_root =
  ("-idas-doc-root", Arg.String (fun d -> idas_doc_root := d), 
   "<dir>  specify the root url for the Sundials IDAS documentation.")
let option_kinsol_doc_root =
  ("-kinsol-doc-root", Arg.String (fun d -> kinsol_doc_root := d), 
   "<dir>  specify the root url for the Sundials KINSOL documentation.")

#if OCAML_3X
let _ =
  let dochtml = new dochtml in
  Odoc_args.add_option option_cvode_doc_root;
  Odoc_args.add_option option_cvodes_doc_root;
  Odoc_args.add_option option_ida_doc_root;
  Odoc_args.add_option option_idas_doc_root;
  Odoc_args.add_option option_kinsol_doc_root;
  Odoc_args.set_doc_generator
    (Some (dochtml :> Odoc_args.doc_generator))
#else
let _ =
  Odoc_args.add_option option_cvode_doc_root;
  Odoc_args.add_option option_cvodes_doc_root;
  Odoc_args.add_option option_ida_doc_root;
  Odoc_args.add_option option_idas_doc_root;
  Odoc_args.add_option option_kinsol_doc_root;
  Odoc_args.extend_html_generator (module Generator : Odoc_gen.Html_functor)
#endif
