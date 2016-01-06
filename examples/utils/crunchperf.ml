(* Separated from perf.ml since updating the formatting shouldn't
   require re-running performance measurements.  *)
let synopsis =
  "crunchperf -c <ocaml> <sundials> <name>

     Combine timing results from previous runs of perf.  <ocaml> should
     list time for OCaml code, <sundials> the time for C code, and <name>
     should identify which example this is.

crunchperf -m <file1> <file2> ...

     Merge multiple log files into one.  Each <file*> should be the
     output of a previous run of crunchperf -c.

crunchperf -s <file>

     Summarize the output of a previous run of crunchperf -m for humans.

crunchperf -S <file>

     Summarize the output of a previous run of crunchperf -m for gnuplot.

crunchperf -r <file> | /bin/sh

     Recover the *.reps files from the output of a previous run of
     crunchperf -m.  Useful when you accidentally made the *.reps files
     out of date but still have the perf.*.log file around.
"

(* Helpers *)

let abbreviate name =
  List.fold_left
    (fun s (pat,rep) -> Str.replace_first (Str.regexp pat) rep s)
    name
    [("^cvode", "cv");
     ("^kinsol", "kin");
     ("^arkode", "ark");
     ("/serial/", "--ser--");
     ("/C_serial/", "--Cser--");
     ("/parallel/", "--par--");
    ]
let expand name =
  List.fold_left
    (fun s (pat,rep) -> Str.replace_first (Str.regexp pat) rep s)
    name
    [("^cv", "cvode");
     ("^kin", "kinsol");
     ("^ark", "arkode");
     ("--ser--", "/serial/");
     ("--Cser--", "/C_serial/");
     ("--par--", "/parallel/");
    ]

let parallel_example s      = Str.string_match (Str.regexp ".*--par--\\|/parallel/.*") s 0
let uses_alternate_solver s = Str.string_match (Str.regexp ".*_alt$") s 0
let uses_nvector_array s    = Str.string_match (Str.regexp ".*_custom$") s 0
let colorof name =
  if parallel_example name then 2
  else if uses_alternate_solver name then 3
  else if uses_nvector_array name then 4
  else 1

let (%) f g x = f (g x)                 (* Function composition *)
let (@@) f x = f x

(* Functions backported from OCaml 4x *)
let mapi f xs =
  let rec go i = function
    | [] -> []
    | x::xs -> f i x :: go (i+1) xs
  in go 0 xs
let iteri f xs =
  let rec go i = function
    | [] -> []
    | x::xs -> f i x; go (i+1) xs
  in go 0 xs

type line = { file : string;
              line : int;
              str : string;             (* Includes trailing '\n'. *)
            }
let get_lines path =
  let file = open_in path in
  let ret = ref [] in
  begin
    try
      while true do
        ret := (input_line file ^ "\n")::!ret
      done
    with End_of_file -> ()
  end;
  mapi (fun i str -> { file = path; line = i+1; str = str })
  @@ List.rev !ret

let comment_line =
  let pat = Str.regexp "[ \t]*\\(#[^\n]*\\)?\n" in
  fun line -> Str.string_match pat line.str 0

let scan_line fmt cont line =
  try Scanf.sscanf line.str fmt cont
  with Scanf.Scan_failure msg ->
    failwith (Printf.sprintf "%s:%d: %s" line.file line.line msg)

type stats = { minimum : float;
               q1 : float;
               median : float;
               q3 : float;
               maximum : float;
               mean : float; }
let analyze dataset =
  let n = Array.length dataset in
  Array.sort compare dataset;
  let maximum = dataset.(n-1)
  and minimum = dataset.(0)
  and q1 = dataset.(n/4)
  and q3 = dataset.(3*n/4)
  and median =
    (* FIXME: take the same kind of care with q1 and q3.  *)
    if (n mod 2) = 0
    then (dataset.(n/2) +. dataset.(n/2-1)) /. 2.
    else dataset.(n/2)
  in
  let total = Array.fold_left (+.) 0. dataset in
  let mean = total /. float_of_int n in
  { mean=mean; minimum=minimum; q1=q1; median=median; q3=q3; maximum=maximum }

(* Input/Output Formats *)
(* OCaml 3.12 allows %.2f in printf but not scanf.  We need two
   essentially identical format strings :-(  *)
let header_with_id = "# ID\treps\tC med.\tOCaml\tC\tOCaml/C\tname\tcategory\n"
let print_fmt_with_id =
  ("%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%d\n"
   : ('a, 'use, 'b, 'c, 'd, 'e) format6)
let scan_fmt_with_id =
  ("%d\t%d\t%f\t%f\t%f\t%f\t%s\t%d\n"
   : ('a, 'use, 'b, 'c, 'd, 'e) format6)

let header_no_id = "# reps\tC med.\tOCaml\tC\tOCaml/C\tname\tcategory\n"
let print_fmt_no_id =
  ("%d\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%d\n"
   : ('a, 'use, 'b, 'c, 'd, 'e) format6)
let scan_fmt_no_id =
  ("%d\t%f\t%f\t%f\t%f\t%s\t%d\n"
   : ('a, 'use, 'b, 'c, 'd, 'e) format6)

let header_for_gnuplot = "# ID\treps\tOCaml\tC\tOCaml/C\tname\tcategory\n"
let print_fmt_for_gnuplot =
  ("%d\t%d\t%.2f\t%.2f\t%.2f\t%s\t%d\n"
   : ('a, 'use, 'b, 'c, 'd, 'e) format6)

let header_for_humans = "# reps\tOCaml\tC\tOCaml/C\tname\n"
let print_fmt_for_humans =
  ("%d\t%.2f\t%.2f\t%.2f\t%s\n"
   : ('a, 'use, 'b, 'c, 'd, 'e) format6)

(* Statistical Analysis *)
type record = { reps : int;
                mutable ml_times : float list;
                mutable c_times : float list;
              }
let load ?(records=Hashtbl.create 10) path =
  let insert _ reps cmed ml c _ name _ =
    try let r = Hashtbl.find records name in
      if reps <> r.reps then failwith (path ^ ": inconsistent reps field");
      r.ml_times <- ml::r.ml_times;
      r.c_times  <- c::r.c_times;
    with Not_found -> Hashtbl.add records name { ml_times = [ml];
                                                 c_times = [c];
                                                 reps = reps; }
  in
  let lines = get_lines path in
  let expect_ids =
    match lines with
    | l::_ -> header_with_id = l.str
    | _ -> false
  in
  List.iter (fun line ->
      if expect_ids
      then scan_line scan_fmt_with_id insert line
      else scan_line scan_fmt_no_id (insert (-1)) line
    )
    (List.filter (not % comment_line) lines);
  records

let sorted_assocs records =
  let ls = ref [] in
  Hashtbl.iter (fun name record -> ls := (name, record)::!ls) records;
  List.sort compare !ls


(* Main routines *)

let combine ocaml sundials name =
  let load path =
    let lines = get_lines path in
    let reps =
      try scan_line "# NUM_REPS = %d\n" (fun r -> r) (List.hd lines)
      with End_of_file | Failure "hd" ->
        failwith ("Input file " ^ path ^ " contains no data")
    in
    let times = ref [] in
    List.iter (scan_line "%f\n" (fun f -> times := f::!times))
    @@ List.filter (not % comment_line) lines;
    if !times = [] then failwith ("Input file " ^ path ^ " contains no data");
    reps, Array.of_list !times
  in
  let ml_reps, ml_times = load ocaml in
  let c_reps, c_times = load sundials in
  if ml_reps <> c_reps then
    Printf.fprintf stderr "Warning: NUM_REPS don't match in %s and %s"
      sundials ocaml;
  if Array.length c_times <> Array.length ml_times then
    failwith (Printf.sprintf "different numbers of data points in %s and %s"
                ocaml sundials);
  print_string header_no_id;
  let c_median = (analyze c_times).median in
  for i = 0 to min (Array.length c_times) (Array.length ml_times) - 1 do
    Printf.printf print_fmt_no_id
      c_reps c_median ml_times.(i) c_times.(i) (ml_times.(i) /. c_times.(i))
      (abbreviate name) (colorof name)
  done

let merge paths =
  let records = Hashtbl.create 100 in
  List.iter (fun p -> ignore (load ~records:records p)) paths;
  print_string header_with_id;
  let records = Array.of_list (sorted_assocs records) in
  let n = Array.length records in
  for id = 0 to n - 1 do
    let name, record = records.(id) in
    let ml = Array.of_list record.ml_times in
    let c = Array.of_list record.c_times in
    let c_median  = (analyze c).median in
    for j = 0 to Array.length c - 1 do
      Printf.printf print_fmt_with_id
        id record.reps c_median ml.(j) c.(j) (ml.(j) /. c.(j))
        name (colorof name)
    done;
    if id < n-1 then print_string "\n\n"
  done

let summarize gnuplot path =
  let assocs = sorted_assocs (load path) in
  Printf.printf "# All numbers except reps show median values.\n";
  if gnuplot
  then print_string header_with_id
  else print_string header_for_humans;
  let _ =
    iteri (fun id (name, record) ->
        let median ls = (analyze (Array.of_list ls)).median in
        let c  = median record.c_times in
        let ml = median record.ml_times in
        let ratio = median (List.map2 (/.) record.ml_times record.c_times) in
        if gnuplot
        then Printf.printf print_fmt_for_gnuplot
               id record.reps ml c ratio name (colorof name)
        else Printf.printf print_fmt_for_humans
               record.reps ml c ratio name)
      assocs
  in
  ()

let recover_reps file =
  let isdir s =
    try Sys.is_directory s
    with _ -> false
  in
  List.iter (fun modname ->
      List.iter (fun nvtype ->
          if not (isdir modname && isdir (Filename.concat modname nvtype))
          then failwith "You don't seem to be in the examples/ directory."
        )
        ["serial"; "parallel"]
    )
    ["cvode"; "cvodes"; "ida"; "idas"; "kinsol"]
  ;
  List.iter
    (fun (name, record) ->
       Printf.printf "echo \"NUM_REPS=%d\" > %s.reps;\n"
         record.reps
         (expand name)
    )
  @@ sorted_assocs @@ load file

let _ =
  match Array.to_list Sys.argv with
  | [_;"-c";ocaml;sundials;name] ->
    combine ocaml sundials name
  | [_;"-s";file] ->
    summarize false file
  | [_;"-S";file] ->
    summarize true file
  | _::"-m"::files -> merge files
  | [_;"-r";file] ->
    recover_reps file
  | _ ->
    print_string synopsis;
    exit 0

