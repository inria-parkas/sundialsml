(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2021 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

[@@@ocaml.warning "-7"]

(**
 Custom tags for the ocamldoc comments:
    @sundials       link to Sundials documentation
    @profiler       link to Sundials Profiler module documentation
    @logger         link to Sundials Logger module documentation
    @context        link to Sundials Context module documentation
    @cvode          link to Sundials CVODE documentation
    @cvodes         link to Sundials CVODES documentation
    @cvodes_quad    link to Sundials CVODES Quadrature documentation
    @cvodes_sens    link to Sundials CVODES Sensitivity documentation
    @cvodes_adj     link to Sundials CVODES Adjoint documentation
    @arkode         link to Sundials ARKODE documentation
    @arkode_ark     link to Sundials ARKODE ARKStep documentation
    @arkode_erk     link to Sundials ARKODE ERKStep documentation
    @arkode_sprk    link to Sundials ARKODE SPRKStep documentation
    @arkode_mri     link to Sundials ARKODE MRIStep documentation
    @arkode_bt      link to Sundials ARKODE Butcher Table documentation
    @arkode_user    link to Sundials ARKODE User-supplied functions documentation
    @arkode_precond link to Sundials ARKODE Preconditioners documentation
    @arkode_innerstepper link to Sundials ARKODE Innerstepper documentation
    @arkode_coupling link to Sundials ARKODE Coupling documentation
    @ida            link to Sundials IDA documentation
    @idas           link to Sundials IDAS documentation
    @idas_quad      link to Sundials IDAS Quadrature documentation
    @idas_sens      link to Sundials IDAS Sensitivity documentation
    @idas_adj       link to Sundials IDAS Adjoint documentation
    @kinsol         link to Sundials KINSOL documentation
    @nvector        link to Sundials NVECTOR documentation
    @matrix         link to Sundials MATRIX documentation
    @matrix_data    link to Sundials MATRIX data-specific documentation
    @linsol         link to Sundials Linear Solver documentation
    @linsol_module  link to Sundials Linear Solver module documentation
    @nonlinsol         link to Sundials Nonlinear Solver documentation
    @nonlinsol_module  link to Sundials Nonlinear Solver module documentation
    @nodoc          No link is possible
 *)

let sundials_doc_root = ref SUNDIALS_DOC_ROOT

let mathjax_url = ref MATHJAX_URL (* directory containing MathJax.js *)

let bp = Printf.bprintf
let bs = Buffer.add_string

type custom_type =
    Simple of (string -> string)
  | Full of (Buffer.t -> Odoc_info.text -> unit)

let sundials_link div_class doc_root page anchor title =
  Printf.sprintf
    "<li class=\"sundials %s\">\
      <span class=\"seesundials\">See Sundials: </span>\
      <a href=\"%s%s%s\">%s</a></li>"
    div_class doc_root page anchor title

module Generator (G : Odoc_html.Html_generator) =
struct
  class html =
  object(self)
    inherit G.html as super

    val rex = Str.regexp "<\\([^#>]*\\)\\(#[^)]*\\)?> \\(.*\\)"

    val variables = [
      ("version", let major, minor, patch, _binding = Sundials_configuration.version in
                  Printf.sprintf "%d.%d.%d" major minor patch)
    ]

    method private split_text ?(page="Usage/index.html") (t:Odoc_info.text) =
      let s = Odoc_info.text_string_of_text t in
      if not (Str.string_match rex s 0) then
        let first_word = List.hd (Str.split (Str.regexp " +") s) in
        (page, "#c." ^ first_word, s)
      else
        let page = Str.matched_group 1 s
        and anchor = try
              Str.matched_group 2 s
            with Not_found -> ""
        and title = Str.matched_group 3 s
        in
        (page, anchor, title)

    method private html_of_missing t =
      let (_page, _anchor, title) = self#split_text t in
      Printf.sprintf
        "<li class=\"sundials nodoc\">\
          <span class=\"seesundials\">See Sundials: </span>%s\
         </li>"
        title

    method private html_of_sundials t =
      let (page, anchor, title) = self#split_text t in
      sundials_link "sundials" (!sundials_doc_root ^ "sundials/") page anchor title

    method private html_of_profiler t =
      let (page, anchor, title) = self#split_text ~page:"Profiling_link.html" t in
      sundials_link "sundials" (!sundials_doc_root ^ "sundials/") page anchor title

    method private html_of_logger t =
      let (page, anchor, title) = self#split_text ~page:"Logging_link.html" t in
      sundials_link "sundials" (!sundials_doc_root ^ "sundials/") page anchor title

    method private html_of_context t =
      let (page, anchor, title) = self#split_text ~page:"SUNContext_link.html" t in
      sundials_link "sundials" (!sundials_doc_root ^ "sundials/") page anchor title

    method private html_of_cvode t =
      let (page, anchor, title) = self#split_text t in
      sundials_link "cvode" (!sundials_doc_root ^ "cvode/") page anchor title

    method private html_of_cvodes t =
      let (page, anchor, title) = self#split_text ~page:"Usage/SIM.html" t in
      sundials_link "cvodes" (!sundials_doc_root ^ "cvodes/") page anchor title

    method private html_of_cvodes_quad t =
      let (page, anchor, title) = self#split_text ~page:"Usage/SIM.html" t in
      sundials_link "cvodes" (!sundials_doc_root ^ "cvodes/") page anchor title

    method private html_of_cvodes_sens t =
      let (page, anchor, title) = self#split_text ~page:"Usage/FSA.html" t in
      sundials_link "cvodes" (!sundials_doc_root ^ "cvodes/") page anchor title

    method private html_of_cvodes_adj t =
      let (page, anchor, title) = self#split_text ~page:"Usage/ADJ.html" t in
      sundials_link "cvodes" (!sundials_doc_root ^ "cvodes/") page anchor title

    method private html_of_arkode t =
      let (page, anchor, title) = self#split_text t in
      sundials_link "arkode" (!sundials_doc_root ^ "arkode/") page anchor title

    method private html_of_arkode_ark t =
      let (page, anchor, title) =
        self#split_text ~page:"Usage/ARKStep_c_interface/User_callable.html" t
      in
      sundials_link "arkode" (!sundials_doc_root ^ "arkode/") page anchor title

    method private html_of_arkode_erk t =
      let (page, anchor, title) =
        self#split_text ~page:"Usage/ERKStep_c_interface/User_callable.html" t
      in
      sundials_link "arkode" (!sundials_doc_root ^ "arkode/") page anchor title

    method private html_of_arkode_sprk t =
      let (page, anchor, title) =
        self#split_text ~page:"Usage/SPRKStep_c_interface/User_callable.html" t
      in
      sundials_link "arkode" (!sundials_doc_root ^ "arkode/") page anchor title

    method private html_of_arkode_mri t =
      let (page, anchor, title) =
        self#split_text ~page:"Usage/MRIStep_c_interface/User_callable.html" t
      in
      sundials_link "arkode" (!sundials_doc_root ^ "arkode/") page anchor title

    method private html_of_arkode_bt t =
      let (page, anchor, title) =
        self#split_text ~page:"ARKodeButcherTable_link.html" t
      in
      sundials_link "arkode" (!sundials_doc_root ^ "arkode/") page anchor title

    method private html_of_arkode_user t =
      let (page, anchor, title) =
        self#split_text ~page:"Usage/User_supplied.html" t
      in
      sundials_link "arkode" (!sundials_doc_root ^ "arkode/") page anchor title

    method private html_of_arkode_precond t =
      let (page, anchor, title) =
        self#split_text ~page:"Usage/ARKStep_c_interface/Preconditioners.html" t
      in
      sundials_link "arkode" (!sundials_doc_root ^ "arkode/") page anchor title

    method private html_of_arkode_innerstepper t =
      let (page, anchor, title) =
        self#split_text ~page:"Usage/MRIStep_c_interface/Custom_Inner_Stepper/Description.html" t
      in
      sundials_link "arkode" (!sundials_doc_root ^ "arkode/") page anchor title

    method private html_of_arkode_coupling t =
      let (page, anchor, title) =
        self#split_text ~page:"Usage/MRIStep_c_interface/MRIStepCoupling.html" t
      in
      sundials_link "arkode" (!sundials_doc_root ^ "arkode/") page anchor title

    method private html_of_ida t =
      let (page, anchor, title) = self#split_text t in
      sundials_link "ida" (!sundials_doc_root ^ "ida/") page anchor title

    method private html_of_idas t =
      let (page, anchor, title) = self#split_text ~page:"Usage/SIM.html" t in
      sundials_link "idas" (!sundials_doc_root ^ "idas/") page anchor title

    method private html_of_idas_quad t =
      let (page, anchor, title) = self#split_text ~page:"Usage/SIM.html" t in
      sundials_link "idas" (!sundials_doc_root ^ "idas/") page anchor title

    method private html_of_idas_sens t =
      let (page, anchor, title) = self#split_text ~page:"Usage/FSA.html" t in
      sundials_link "idas" (!sundials_doc_root ^ "idas/") page anchor title

    method private html_of_idas_adj t =
      let (page, anchor, title) = self#split_text ~page:"Usage/ADJ.html" t in
      sundials_link "idas" (!sundials_doc_root ^ "idas/") page anchor title

    method private html_of_kinsol t =
      let (page, anchor, title) = self#split_text t in
      sundials_link "kinsol" (!sundials_doc_root ^ "kinsol/") page anchor title

    method private html_of_nvector t =
      let (page, anchor, title) = self#split_text ~page:"NVector_links.html" t in
      sundials_link "nvector" (!sundials_doc_root ^ "nvectors/") page anchor title

    method private html_of_matrix t =
      let (page, anchor, title) = self#split_text ~page:"SUNMatrix_API_link.html" t in
      sundials_link "matrix" (!sundials_doc_root ^ "sunmatrix/") page anchor title

    method private html_of_matrix_data t =
      let (page, anchor, title) = self#split_text ~page:"SUNMatrix_links.html" t in
      sundials_link "matrix" (!sundials_doc_root ^ "sunmatrix/") page anchor title

    method private html_of_linsol t =
      let (page, anchor, title) = self#split_text ~page:"SUNLinSol_API_link.html" t in
      sundials_link "linsol" (!sundials_doc_root ^ "sunlinsol/") page anchor title

    method private html_of_linsol_module t =
      let (page, anchor, title) = self#split_text ~page:"SUNLinSol_links.html" t in
      sundials_link "linsol" (!sundials_doc_root ^ "sunlinsol/") page anchor title

    method private html_of_nonlinsol t =
      let (page, anchor, title) = self#split_text ~page:"SUNNonlinSol_API_link.html" t in
      sundials_link "nonlinsol" (!sundials_doc_root ^ "sunnonlinsol/") page anchor title

    method private html_of_nonlinsol_module t =
      let (page, anchor, title) = self#split_text ~page:"SUNNonlinSol_links.html" t in
      sundials_link "nonlinsol" (!sundials_doc_root ^ "sunnonlinsol/") page anchor title

    val divrex = Str.regexp " *\\(open\\|close\\) *\\(.*\\)"

    method private html_of_div s =
      let baddiv s = (Odoc_info.warning (Printf.sprintf
            "div must be followed by 'open' or 'close', not '%s'!" s); "") in
      if not (Str.string_match divrex s 0) then baddiv s
      else
        match Str.matched_group 1 s with
        | "open" ->
            let attrs = try Str.matched_group 2 s with Not_found -> "" in
            Printf.sprintf "<div %s>" attrs
        | "close" -> "</div>"
        | s -> baddiv s

    method private html_of_var s =
      let var = Str.replace_first (Str.regexp " +$") ""
                  (Str.replace_first (Str.regexp "^ +") "" s) in
      try
        List.assoc var variables
      with Not_found ->
        (Odoc_info.warning (Printf.sprintf "Variable '%s' is not defined." var); "")

    method private html_of_img s =
      Printf.sprintf "<a href=\"%s\"><img src=\"%s\"></a>" s s

    method private html_of_openfile s =
      let var = Str.replace_first (Str.regexp " +$") ""
                  (Str.replace_first (Str.regexp "^ +") "" s) in
      Printf.sprintf "<a href=\"%s\">%s</a>" var var

    method private html_of_cconst s =
      let var = Str.replace_first (Str.regexp " +$") ""
                  (Str.replace_first (Str.regexp "^ +") "" s) in
      Printf.sprintf "<span class=\"cconst\">(%s)</span>" var

    method private html_of_color s =
      let ss = Str.bounded_split (Str.regexp "[ \t\n]+") s 2 in
      match ss with
      | [x] ->
          (Odoc_info.warning (Printf.sprintf "No color given ('%s')." x); x)
      | [color; text] ->
          Printf.sprintf "<span style=\"color: %s;\">%s</span>" color text
      | _ -> assert false

    method html_of_raised_exceptions b l =
      match l with
        [] -> ()
      | (s, t) :: [] ->
          bs b "<div class=\"raisedexceptions\">";
          bp b "<span class=\"raises\">%s</span> <code>%s</code> "
            Odoc_messages.raises
            s;
          self#html_of_text b t;
          bs b "</div>\n"
      | _ ->
          bs b "<div class=\"raisedexceptions\">";
          bp b "<span class=\"raises\">%s</span><ul>" Odoc_messages.raises;
          List.iter
            (fun (ex, desc) ->
              bp b "<li><code>%s</code> " ex ;
              self#html_of_text b desc;
              bs b "</li>\n"
            )
            l;
          bs b "</ul></div>\n"

    method html_of_author_list b l =
      match l with
        [] -> ()
      | _ ->
          bp b "<div class=\"authors\">";
          bp b "<b>%s:</b> " Odoc_messages.authors;
          self#html_of_text b [Odoc_info.Raw (String.concat ", " l)];
          bs b "</div>\n"

    val mutable custom_functions =
      ([] : (string * custom_type) list)

    method private html_of_warning b t =
      bs b "<div class=\"warningbox\">";
      self#html_of_text b t;
      bs b "</div>"

    method private html_of_custom_text b tag text =
      try
        match List.assoc tag custom_functions, text with
        | (Simple f, [Odoc_info.Raw s]) -> Buffer.add_string b (f s)
        | (Simple _, _) ->
            Odoc_info.warning (Printf.sprintf 
              "custom tags (%s) must be followed by plain text." tag)
        | (Full f, _) -> f b text
      with
        Not_found -> Odoc_info.warning (Odoc_messages.tag_not_handled tag)

    (* Import MathJax (http://www.mathjax.org/) to render mathematics in
       function comments. *)
    method init_style =
      super#init_style;
      style <- style ^
                "<script type=\"text/x-mathjax-config\">\n" ^
                "   MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$']]}});\n" ^
                "</script>" ^
                "<script type=\"text/javascript\"\n" ^
                Printf.sprintf
                  "        src=\"%s/MathJax.js?config=TeX-AMS-MML_HTMLorMML\">\n"
                  !mathjax_url ^
                "</script>\n"

    method html_of_Latex b s =
      Buffer.add_string b s

    initializer
      tag_functions <- ("sundials", self#html_of_sundials) :: tag_functions;
      tag_functions <- ("profiler", self#html_of_profiler) :: tag_functions;
      tag_functions <- ("logger", self#html_of_logger) :: tag_functions;
      tag_functions <- ("context", self#html_of_context) :: tag_functions;
      tag_functions <- ("cvode",    self#html_of_cvode) :: tag_functions;
      tag_functions <- ("cvodes",   self#html_of_cvodes) :: tag_functions;
      tag_functions <- ("cvodes_quad", self#html_of_cvodes_quad) :: tag_functions;
      tag_functions <- ("cvodes_sens", self#html_of_cvodes_sens) :: tag_functions;
      tag_functions <- ("cvodes_adj",  self#html_of_cvodes_adj)  :: tag_functions;
      tag_functions <- ("arkode",   self#html_of_arkode) :: tag_functions;
      tag_functions <- ("arkode_ark", self#html_of_arkode_ark) :: tag_functions;
      tag_functions <- ("arkode_erk", self#html_of_arkode_erk) :: tag_functions;
      tag_functions <- ("arkode_sprk", self#html_of_arkode_sprk) :: tag_functions;
      tag_functions <- ("arkode_mri", self#html_of_arkode_mri) :: tag_functions;
      tag_functions <- ("arkode_bt",  self#html_of_arkode_bt) :: tag_functions;
      tag_functions <- ("arkode_user",self#html_of_arkode_user) :: tag_functions;
      tag_functions <- ("arkode_precond",self#html_of_arkode_precond) :: tag_functions;
      tag_functions <- ("arkode_innerstepper",self#html_of_arkode_innerstepper) :: tag_functions;
      tag_functions <- ("arkode_coupling",self#html_of_arkode_coupling) :: tag_functions;
      tag_functions <- ("ida",      self#html_of_ida) :: tag_functions;
      tag_functions <- ("idas",     self#html_of_idas) :: tag_functions;
      tag_functions <- ("idas_quad", self#html_of_idas_quad) :: tag_functions;
      tag_functions <- ("idas_sens", self#html_of_idas_sens) :: tag_functions;
      tag_functions <- ("idas_adj",  self#html_of_idas_adj)  :: tag_functions;
      tag_functions <- ("kinsol",   self#html_of_kinsol) :: tag_functions;
      tag_functions <- ("nvector",  self#html_of_nvector) :: tag_functions;
      tag_functions <- ("matrix",   self#html_of_matrix) :: tag_functions;
      tag_functions <- ("matrix_data", self#html_of_matrix_data) :: tag_functions;
      tag_functions <- ("linsol",        self#html_of_linsol) :: tag_functions;
      tag_functions <- ("linsol_module", self#html_of_linsol_module) :: tag_functions;
      tag_functions <- ("nonlinsol",        self#html_of_nonlinsol) :: tag_functions;
      tag_functions <- ("nonlinsol_module", self#html_of_nonlinsol_module) :: tag_functions;
      tag_functions <- ("nodoc",   self#html_of_missing) :: tag_functions;

      custom_functions <- ("div",      Simple self#html_of_div)      ::
                          ("var",      Simple self#html_of_var)      ::
                          ("color",    Simple self#html_of_color)    ::
                          ("img",      Simple self#html_of_img)      ::
                          ("cconst",   Simple self#html_of_cconst)   ::
                          ("openfile", Simple self#html_of_openfile) ::
                          ("warning", Full self#html_of_warning)     ::
                          custom_functions

  end
end

let _  = Odoc_html.charset := "utf-8"

let option_sundials_doc_root =
  ("-sundials-doc-root", Arg.String (fun d -> sundials_doc_root := d), 
   "<dir>  specify the root url for the Sundials documentation.")
let option_mathjax_url =
  ("-mathjax", Arg.String (fun d -> mathjax_url := d), 
   "<url>  specify the root url for MathJax.")

let _ =
  Odoc_args.add_option option_sundials_doc_root;
  Odoc_args.add_option option_mathjax_url;
  Odoc_args.extend_html_generator (module Generator : Odoc_gen.Html_functor)

