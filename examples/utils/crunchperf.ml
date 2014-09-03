(* Separated from perf.ml since updating the formatting shouldn't
   require re-running performance measurements.  *)
let synopsis =
"crunchperf -c <ocaml> <sundials> <name>

     Combine two performance measurements from previous runs of
     perf.  <ocaml> should list the OCaml code's run time, <sundials>
     the C code's, and <name> should identify which example this is.

crunchperf -s <file>

     Summarize the output of a previous run of crunchperf -p for humans.

crunchperf -S <file>

     Summarize the output of a previous run of crunchperf -p for gnuplot.
"

let abbreviate name =
  List.fold_left
    (fun s (pat,rep) -> Str.replace_first (Str.regexp pat) rep s)
    name
    [("^cvode", "cv");
     ("^kinsol", "kin");
     ("/serial/", "/ser/");
     ("/parallel/", "/par/");
    ]
let expand name =
  List.fold_left
    (fun s (pat,rep) -> Str.replace_first (Str.regexp pat) rep s)
    name
    [("^cv", "cvode");
     ("^kin", "kinsol");
     ("/ser/", "/serial/");
     ("/par/", "/parallel/");
    ]

(* Helpers *)
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

(* Discard a comment or empty line.  We can't use input_line, as it
   interferes with fscanf's buffering.  *)
let discard_hash_line path file msg_when_fail =
  let rec read_till_newline file =
    if Scanf.fscanf file "%c" (fun c -> c) <> '\n'
    then read_till_newline file
  in
  try Scanf.fscanf file " #" (); read_till_newline file
  with End_of_file -> ()
     | Scanf.Scan_failure _ ->
       try Scanf.fscanf file " \n" ()
       with End_of_file -> ()
          | Scanf.Scan_failure _ ->
            Scanf.fscanf file "%l" (fun line ->
                failwith (Printf.sprintf "%s:%d: %s" path line msg_when_fail))

let rec scan_line path file fmt cont =
  try Scanf.fscanf file fmt cont
  with Scanf.Scan_failure msg ->
    discard_hash_line path file msg;
    scan_line path file fmt cont

let combine ocaml sundials name =
  let load path =
    let file = open_in path in
    let reps =
      try Scanf.fscanf file "# NUM_REPS = %d\n" (fun r -> r)
      with End_of_file -> failwith ("Input file " ^ path ^ " contains no data")
    in
    let times = ref [] in
    begin
      try while true do
            scan_line path file "%.2f\n" (fun f -> times := f::!times)
          done;
      with End_of_file -> ()
    end;
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
  Printf.printf "# reps\tC med.\tOCaml\tC\tOCaml/C\tname\n";
  let c_median = (analyze c_times).median in
  for i = 0 to Array.length c_times - 1 do
    Printf.printf "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n"
      c_reps c_median ml_times.(i) c_times.(i) (ml_times.(i) /. c_times.(i))
      (abbreviate name)
  done

type record = { reps : int;
                mutable ml_times : float list;
                mutable c_times : float list;
              }
let summarize abbrev path =
  let file = open_in path in
  let records = Hashtbl.create 10 in
  let insert reps _ ml c _ name =
    try let r = Hashtbl.find records name in
      if reps <> r.reps then failwith (path ^ ": inconsistent reps field");
      r.ml_times <- ml::r.ml_times;
      r.c_times  <- c::r.c_times;
    with Not_found -> Hashtbl.add records name { ml_times = [ml];
                                                 c_times = [c];
                                                 reps = reps; }
  in
  begin
    try while true do
          scan_line path file "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n" insert
        done;
    with End_of_file -> ()
  end;
  Printf.printf "# All numbers except reps show median values.\n";
  Printf.printf "# reps\tOCaml\tC\tOCaml/C\tname\n";
  let assocs = ref [] in
  Hashtbl.iter (fun name record -> assocs := (name, record)::!assocs) records;
  let assocs = List.sort compare !assocs in
  List.iter (fun (name, record) ->
      let median ls = (analyze (Array.of_list ls)).median in
      let c  = median record.c_times in
      let ml = median record.ml_times in
      let ratio = median (List.map2 (/.) record.ml_times record.c_times) in
      Printf.printf "%d\t%.2f\t%.2f\t%.2f\t%s\n"
        record.reps ml c ratio (if abbrev then name else expand name))
    assocs

let _ =
  match Array.to_list Sys.argv with
  | [_;"-c";ocaml;sundials;name] ->
    combine ocaml sundials name
  | [_;"-s";file] ->
    summarize false file
  | [_;"-S";file] ->
    summarize true file
  | _ ->
    print_string synopsis;
    exit 0

