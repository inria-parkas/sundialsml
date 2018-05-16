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

crunchperf -s [-u [<confidence>]] <file>

     Summarize the output of a previous run of crunchperf -m for humans.
     If given -u, compute a range of OCaml/C ratios that are statistically
     significant (requires GNU Octave; <confidence> defaults to 0.005).

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
     ("/C_serial/", "--C_ser--");
     ("/parallel/", "--par--");
     ("/C_parallel/", "--C_par--");
     ("/C_openmp/", "--C_omp--");
    ]
let expand name =
  List.fold_left
    (fun s (pat,rep) -> Str.replace_first (Str.regexp pat) rep s)
    name
    [("^cv", "cvode");
     ("^kin", "kinsol");
     ("^ark", "arkode");
     ("--ser--", "/serial/");
     ("--C_ser--", "/C_serial/");
     ("--par--", "/parallel/");
     ("--C_par--", "/C_parallel/");
     ("--C_omp--", "/C_openmp/");
    ]

let parallel_example s      = Str.string_match (Str.regexp ".*--\\(C_\\)?par--\\|/parallel/.*") s 0
let uses_alternate_solver s = Str.string_match (Str.regexp ".*_alt$") s 0
let uses_nvector_array s    = Str.string_match (Str.regexp ".*_custom$") s 0
let colorof name =
  if parallel_example name then 2
  else if uses_alternate_solver name then 3
  else if uses_nvector_array name then 4
  else 1

let (%) f g x = f (g x)                 (* Function composition *)
let (@@) f x = f x

let string_of_float_array ?(sep=";") a =
  "[" ^ String.concat sep (Array.to_list @@ Array.map string_of_float a) ^ "]"

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

let array_iteri f xs =
  for i = 0 to Array.length xs - 1 do
    f i xs.(i)
  done

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

let header_for_humans_with_range = "# reps\tOCaml med.\tC med.\tOCaml/C lower\tOCaml/C upper\tname\n"
let print_fmt_for_humans_with_range =
  ("%d\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n"
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

let rank_sums xs ys =
  let xs = Array.map (fun x -> (x, ref nan)) xs
  and ys = Array.map (fun y -> (y, ref nan)) ys
  in
  let samples = Array.append xs ys in
  Array.sort compare samples;
  let rec assign_ranks i =
    if i < Array.length samples
    then
      let rec probe j =
        if j < Array.length samples && fst samples.(i) = fst samples.(j)
        then probe (j+1)
        else j
      in
      let next_i = probe (i+1) in
      let mean_rank = float_of_int ((i+1) + next_i) /. 2. in
      for k = i to next_i - 1 do
        snd samples.(k) := mean_rank
      done;
      assign_ranks next_i
  in
  assign_ranks 0;
  let rank_sum = Array.fold_left (fun a (_, r) -> a +. !r) 0. in
  (rank_sum xs, rank_sum ys)

(* FIXME: This is really slow because it re-launches octave over and
   over.  *)
let u_test ?(side="<>") xs ys =
  (* Returns p-value of Mann-Whitney test with the null hypothesis
     P(X > Y) = P(Y < X).  *)
  let octave_path =
    try Unix.getenv "OCTAVE_CLI"
    with Not_found -> "octave-cli"
  and script =
    Printf.sprintf
      "[p,_] = u_test(%s, %s, \"%s\"); disp(p)"
      (string_of_float_array ~sep:"," xs)
      (string_of_float_array ~sep:"," ys)
      side
  in
  let command = Printf.sprintf "%s --eval '%s'" octave_path script in
  let pipe = Unix.open_process_in command
  in
  try
    let ret = Scanf.bscanf (Scanf.Scanning.from_channel pipe)
                           " %f\n" (fun f -> f) in
    close_in pipe;
    ret
  with e -> (close_in pipe; raise e)

let u_interval ?(confidence=0.005) xs ys =
  (* Returns the range of gamma such that, when xs is scaled by gamma,
     the U-test does NOT detect a difference between P(X > Y) and P(Y < X)
     at the given confidence level.  *)
  assert(Array.length xs > 0);
  assert(Array.length ys > 0);
  let xmin = Array.fold_left min xs.(0) xs
  and xmax = Array.fold_left max xs.(0) xs
  and ymin = Array.fold_left min ys.(0) ys
  and ymax = Array.fold_left max ys.(0) ys
  in
  let scale gamma xs = Array.map (fun x -> x *. gamma) xs in
  let p_value gamma = u_test (scale gamma xs) ys in
  let rec bsearch eps pred l h =
    (* Invariant: pred l = false, pred h = true, and pred is
       monotonic.  Return the least number for which pred is true.
       Works on inverted ranges too.  *)
    if abs_float (h -. l) <= eps then h
    else
      let m = (l +. h) /. 2. in
      if pred m then bsearch eps pred l m
      else bsearch eps pred m h
  in
  (* P-value is minimal when one sample dominates the other.  Make
     sure that minimum goes below the desired confidence level.  *)
  let ones  = Array.make (Array.length xs) 1.0
  and zeros = Array.make (Array.length xs) 0.0
  in
  (if u_test ones zeros > confidence
   then raise (Invalid_argument
                 (Printf.sprintf
                    "not enough data points (%d, %d) for confidence level %f"
                    (Array.length xs) (Array.length ys) confidence))
  );
  (if u_test zeros ones > confidence
   then raise (Invalid_argument
                 (Printf.sprintf
                    "not enough data points (%d, %d) for confidence level %f"
                    (Array.length xs) (Array.length ys) confidence))
  );
  let eps = 1.0 +. 1e-10 in
  let gmin = (ymin /. xmax) /. eps
  and gmax = (ymax /. xmin) *. eps in
  (* We need at least one point where statistical significance fails.
     The most likely candidate is where the data crosses over from ys
     being dominant to scaled xs being dominant.  *)
  let rank_sum_gt g =
    let (xr, yr) = rank_sums (scale g xs) ys in
    xr > yr
  in
  assert (not (rank_sum_gt gmin));
  assert (rank_sum_gt gmax);
  let mid = bsearch 1e-3 rank_sum_gt gmin gmax in
  (if p_value mid <= confidence
   then raise (Invalid_argument
                 (Printf.sprintf
                    "p-value always below threshold; not enough data points?"))
  );
  let tol = (gmin +. gmax) *. 1e-3 in
  (bsearch tol (fun g -> p_value g > confidence) gmin mid,
   bsearch tol (fun g -> p_value g <= confidence) mid gmax)

(* Main routines *)

let combine ocaml sundials name =
  let load path =
    let lines = get_lines path in
    let reps =
      try scan_line "# NUM_REPS = %d\n" (fun r -> r) (List.hd lines)
      with End_of_file | Failure _ ->
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

let merge_u paths =
  let records = Hashtbl.create 100 in
  List.iter (fun p -> ignore (load ~records:records p)) paths;
  print_string 

let summarize gnuplot ?confidence path =
  let assocs = sorted_assocs (load path) in
  (match confidence with
   | Some conf ->
      (print_string "# OCaml and C fields show median values.\n";
       Printf.printf "# confidence level = %f\n" conf;
       print_string (if gnuplot then raise (Invalid_argument "not implemented") (* header_with_id_and_range *)
                     else header_for_humans_with_range))
   | None ->
      (print_string "# OCaml, C, and OCaml/C fields show median values.\n";
       print_string (if gnuplot then header_with_id
                     else header_for_humans)));
  let summarize_record id (name, record) =
    let median ls = (analyze ls).median in
    let c = Array.of_list record.c_times in
    let ml = Array.of_list record.ml_times in
    let c_med  = median c in
    let ml_med = median ml in
    let ratio = median (Array.map2 (/.) ml c) in
    match confidence with
    | Some conf ->
       let (low, high) = u_interval ~confidence:conf c ml in
       if gnuplot
       then raise (Invalid_argument "not implemented")
       else Printf.printf print_fmt_for_humans_with_range
                          record.reps ml_med c_med low high name
    | None ->
       if gnuplot
       then Printf.printf print_fmt_for_gnuplot
                          id record.reps ml_med c_med ratio name (colorof name)
       else Printf.printf print_fmt_for_humans
                          record.reps ml_med c_med ratio name
  in
  ignore (iteri summarize_record assocs)

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
  if not !Sys.interactive
  then
    match Array.to_list Sys.argv with
    | [_;"-c";ocaml;sundials;name] ->
       combine ocaml sundials name
    | _::"-s"::rest ->
       (match rest with
        | ["-u";conf;file] ->
           summarize ~confidence:(float_of_string conf) false file
        | ["-u";file] ->
           summarize ~confidence:0.005 false file
        | [file] ->
           summarize false file
        | _ -> print_string synopsis; exit 0)
    | [_;"-S";file] ->
       summarize true file
    | _::"-m"::files -> merge files
    | [_;"-r";file] ->
       recover_reps file
    | _ ->
       print_string synopsis;
       exit 0

