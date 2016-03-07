let synopsis =
"perf -r <min time> <command>

     Determine the number of reps needed for <command> to take at least
     about <min time> (wall-clock time, in seconds; fractional value
     accepted).  The <command> should take time roughly proportional to
     the (integral) numerical value of the environment variable NUM_REPS.

perf -m <reps file> <n> <command>

     Take <n> measurements of the execution time (wall-clock time) of
     <command> with NUM_REPS set to the value recorded in <reps file>.
     <reps file> should contain the output of a previous run of perf -r.
"
let (@@) f x = f x
let init_stop_watch executable args =
  let dev_null =
    try Unix.openfile "/dev/null" [Unix.O_RDWR; Unix.O_TRUNC] 0
    with Unix.Unix_error (Unix.ENOENT, "open", "/dev/null") ->
      (* Windows *)
      try Unix.openfile "NUL" [Unix.O_RDWR; Unix.O_TRUNC] 0
      with Unix.Unix_error (Unix.ENOENT, "open", "NUL") ->
        Printf.fprintf stderr
          "ERROR: Can't find /dev/null or equivalent on your system.  \
           Giving up.";
        exit 2
  in
  let pipe_read, pipe_write = Unix.pipe () in
  let pipe_read = Unix.in_channel_of_descr pipe_read in

  (* The current environment, with NUM_REPS=* moved to index 0. *)
  let env =
    let nlen = String.length "NUM_REPS=" in
    let not_num_reps s =
      String.length s < nlen || String.sub s 0 nlen <> "NUM_REPS="
    in
    Array.of_list ("NUM_REPS=0"::List.filter not_num_reps
                     (Array.to_list (Unix.environment ())))
  in
  let dump_env env =
    Printf.fprintf stderr "Environment:\n";
    Array.iter (fun s -> Printf.fprintf stderr "%s\n" s) env
  in

  (* Use time(1) if one is installed and accepts -f '%e'.  Otherwise,
     use Unix.gettimeofday.  The latter has more overhead (and
     probably less accurate) but is more portable.  *)
  let spawn_wait file args stdin stdout stderr =
    let pid = Unix.create_process_env file args env stdin stdout stderr in
    let cmd () = env.(0) ^ " " ^ String.concat " " (Array.to_list args) in
    let _, status = Unix.waitpid [] pid in
    match status with
    | Unix.WEXITED 0 -> ()
    | Unix.WEXITED n ->
       dump_env env;
       failwith ("Command " ^ cmd () ^ " exited with nonzero status "
                 ^ string_of_int n)
    | Unix.WSIGNALED n ->
       dump_env env;
       failwith ("Command " ^ cmd ()
                 ^ " killed by signal "
                 ^ string_of_int n)
    | Unix.WSTOPPED n ->
       dump_env env;
       failwith ("Command stopped by signal - execution time "
                 ^ "measurement is compromised.")
  in
  let with_gettimeofday =
    let args = Array.append [|executable|] args in
    fun num_reps ->
      env.(0) <- Printf.sprintf "NUM_REPS=%d" num_reps;
      let start = Unix.gettimeofday () in
      spawn_wait executable args dev_null dev_null Unix.stderr;
      let finish = Unix.gettimeofday () in
      finish -. start
  and with_time_command executable args =
    let args = Array.append [|"time"; "-f"; "%e"; executable|] args in
    fun num_reps ->
      env.(0) <- Printf.sprintf "NUM_REPS=%d" num_reps;
      spawn_wait "time" args dev_null dev_null pipe_write;
      Scanf.bscanf (Scanf.Scanning.from_channel pipe_read) "%f\n" (fun f -> f)
  in
  try ignore (with_time_command "true" [||] 1);
    with_time_command executable args
  with _ ->
    prerr_string
      ("Warning: can't find time(1) that accepts -f, measuring time\n" ^
       "with OCaml.  This may be slightly less accurate than time(1).\n");
    flush stderr;
    with_gettimeofday

let determine_reps min_time executable args =
  if min_time < 0.1 then failwith "min time is too small; it should be >= 0.1";
  let measure = init_stop_watch executable (Array.of_list args) in
  let rec go reps =
    let t = measure reps in
    if t >= min_time then reps
    else
      (* Multiply reps by min_time / t.  When t is small, this value
         is likely to be inaccurate, so we put a floor on t.  We also
         want to bump reps by at least 10% of its current value, to
         ensure quick convergence.  *)
      let next_reps =
        int_of_float (ceil (float_of_int
                              reps *. max 1.1 (min_time /. max 0.1 t)))
      in
      if next_reps <= reps then
        failwith @@
        Printf.sprintf "Bug in perf.ml: reps=%d, next_reps=%d, t=%g\n"
          reps next_reps t;
      go next_reps
  in
  let reps = go 1 in
  let quote s = "'" ^ s ^ "'" in
  Printf.printf "NUM_REPS=%d\n\
                 # min time %g[s]\n\
                 # command %s\n"
    reps min_time
    (String.concat " " (List.map quote (executable::args)))

let measure_performance reps n executable args =
  let measure = init_stop_watch executable (Array.of_list args) in
  Printf.printf "# NUM_REPS=%d\n" reps; flush stdout;
  for i = 1 to n do
    let t = measure reps in
    Printf.printf "%.2f\n" t; flush stdout;
  done

let load_reps_file file =
  try Scanf.bscanf (Scanf.Scanning.open_in file) "NUM_REPS = %d\n" (fun r -> r)
  with e -> failwith @@ Printf.sprintf "Reps file %s broken: %s" file
                     @@ Printexc.to_string e

(* Parse args. *)

let read_pos_int s =
  try let n = int_of_string s in
    if n <= 0 then failwith ""
    else n
  with Failure _ ->
    failwith ("Error: expected positive integer but got " ^ s)
let read_pos_float s =
  try let f = float_of_string s in
    if f <= 0. then failwith ""
    else f
  with Failure _ ->
    failwith ("Error: positive real number but got " ^ s)

let _ =
  match Array.to_list Sys.argv with
  | _::"-r"::min_time::ex::args ->
    determine_reps (read_pos_float min_time) ex args
  | _::"-m"::reps_file::n::ex::args ->
    let n = read_pos_int n in
    let reps = load_reps_file reps_file in
    measure_performance reps n ex args
  | _ ->
    print_string synopsis;
    exit 0

