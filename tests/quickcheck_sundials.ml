(* Generators, shrinkers, and pretty-printers for sundials types that are
   shared between IDA and CVODE, as well as some auxiliary functions for
   compiling and running test cases.  Pretty-printers derived in
   pprint_sundials.ml are wrapped or overridden to fine-tune the output.  *)
open Pprint
open Quickcheck
open Sundials
include Pprint_sundials

(* Generators and shrinkers for different kinds of time values.  Time values
   can occur in several places in a script:

   - the query time of SolveNormal

   - the time of the zero of a root function

   - the stop time

   We have to make sure that these types of events do not happen simultaneously
   because otherwise the order in which IDA detects them is dictated by
   floating point error and is unpredictable.

   Therefore, we:

   - discretize the time values, i.e. make them a multiple of a fixed value
     called discrete_unit

   - offset each type of time value by a sub-discrete_unit value, so that
     e.g. the stop time is always distinct from a query time modulo
     discrete_unit.

   Be careful, because this means that the set of valid times is not
   necessarily closed under addition!  For instance, a valid stop time + a
   valid stop time is not a valid stop time.  This is not a big deal, with one
   exception: we need to do a fair amount of arithmetic on query times, so we
   specially require query_time_offs to be zero.

 *)
let discrete_unit = 1.
(* NB: lots of code depend on the fact that query_time_offs is zero! *)
let query_time_offs = 0.
let root_time_offs = discrete_unit /. 2.
let stop_time_offs = discrete_unit /. 4.

(* This must be smaller than any of the offsets.  *)
let time_epsilon = discrete_unit /. 8.

type sign = Positive | NonNegative | ArbitrarySign
(* Returns (gen, shrink), where shrink's first argument gives a *non-strict*
   lower bound on the allowable time values.  *)
let discrete_float_type ?(sign=NonNegative) offset =
  let gen_rank, shrink_rank =
    match sign with
    | Positive -> gen_pos, shrink_pos
    | NonNegative -> gen_nat, shrink_nat
    | ArbitrarySign -> gen_int, shrink_int
  in
  ((fun t0 () -> float_of_int (gen_rank ()) *. discrete_unit +. offset +. t0),
   (fun t0 f ->
     Fstream.map
       (fun x -> float_of_int x *. discrete_unit +. offset +. t0)
       (shrink_rank (int_of_float ((f -. t0 -. offset) /. discrete_unit)))))

let gen_t0, shrink_t0 = discrete_float_type ~sign:ArbitrarySign 0.
let gen_t0, shrink_t0 = gen_t0 0., shrink_t0 0.

let gen_query_time, shrink_query_time = discrete_float_type query_time_offs
let gen_root_time, shrink_root_time = discrete_float_type root_time_offs
let gen_stop_time, shrink_stop_time = discrete_float_type stop_time_offs

(* Similar arrangement for non-time values.  *)
let gen_discrete_float, shrink_discrete_float = discrete_float_type 0.


let gen_carray gen_elem =
  gen_bigarray1 Bigarray.float64 Bigarray.c_layout gen_elem

let shrink_carray ?(shrink_size=true) shrink_elem =
  shrink_bigarray1 Bigarray.float64 Bigarray.c_layout ~shrink_size shrink_elem

let carray_drop_elem : Carray.t -> int -> Carray.t = bigarray1_drop_elem

let pp_carray, dump_carray, show_carray, display_carray,
  print_carray, prerr_carray =
  printers_of_pp (fun ?(prec=0) fmt xs ->
    if !read_write_invariance
    then pp_parens (prec >= Prec.app)
          (pp_array_like Carray.length Bigarray.Array1.get
           "Carray.of_array [|@[<hov>" "@]|]" pp_float) fmt xs
    else pp_array_like Carray.length Bigarray.Array1.get
           "[<@[<hov>" "@]>]" pp_float fmt xs)


let gen_root_event, shrink_root_event =
  gen_shrink_choice [| Roots.NoRoot; Roots.Rising; Roots.Falling |]

let gen_root_info = gen_array_like Roots.make Roots.set gen_root_event

let pp_root_info, dump_root_info, show_root_info, display_root_info,
  print_root_info, prerr_root_info
    =
  printers_of_pp (fun ?(prec=0) fmt xs ->
    let get a i =
    (* If Roots is used improperly or if there's a bug in the binding, a root
       info array can contain garbage that doesn't correspond to any of NoRoot,
       Rising, or Falling.  Such values trigger a Failure in Roots.get.  *)
      try show_root_event (Roots.get a i)
      with Failure _ -> "<garbage>"
    in
    if !read_write_invariance
    then pp_parens (prec >= Prec.app)
           (pp_array_like Roots.length get "Roots.of_array [|@[<hov>" "@]|]"
              pp_string_noquote)
           fmt xs
    else pp_array_like Roots.length get "[<" ">]" pp_string_noquote fmt xs)

let gen_root_direction, shrink_root_direction =
  gen_shrink_choice [| RootDirs.Increasing; RootDirs.Decreasing;
                       RootDirs.IncreasingOrDecreasing |]

let gen_root_dirs =
  gen_array_like RootDirs.make RootDirs.set gen_root_direction

let shrink_root_dirs =
  shrink_array_like RootDirs.make RootDirs.length RootDirs.get RootDirs.set
    shrink_root_direction

let pp_root_dirs, dump_root_dirs, show_root_dirs, display_root_dirs
  , print_root_dirs, prerr_root_dirs =
  printers_of_pp
    (fun ?(prec=0) fmt xs ->
      let get a i =
        (* If RootDirs is used improperly or if there's a bug in the binding, a
           root info array can contain garbage that doesn't correspond to any
           of Increasing, Decreasing, or IncreasingOrDecreasing.  Such values
           trigger a Failure in RootDirs.get.  *)
        try show_root_direction (RootDirs.get a i)
        with Failure _ -> "<garbage>"
      in
      if !read_write_invariance
      then pp_parens (prec >= Prec.app)
             (pp_array_like RootDirs.length get "RootDirs.of_array [|" "|]"
                pp_string_noquote)
             fmt xs
      else pp_array_like RootDirs.length get "[<" ">]"
             pp_string_noquote fmt xs)

type result =
    Unit
  | Int of int | Float of float
  | Any
  | Type of result
  | Aggr of result list
  | Carray of Carray.t    (* NB: always copy the array! *)
  | SolverResult of solver_result
  | Exn of exn
  | RootInfo of Roots.t
deriving (pretty ~alias:(Carray.t = carray,
                         Roots.t = root_info)
                 ~optional:(Int, Float))

let pp_result, dump_result, show_result, display_result,
  print_result, prerr_result =
  printers_of_pp (fun ?(prec=0) fmt x ->
      if !read_write_invariance
      then pp_result ~prec fmt x (* use the derived one *)
      else
        match x with
        | Unit -> pp_string_noquote fmt "()"
        | Any ->  pp_string_noquote fmt "_"
        | x -> pp_result ~prec fmt x (* use the derived one *)
    )

(* Use this instead of the constructor Carray.  *)
let carray a = Carray (Carray.of_carray a)

(* Check if r1 is a valid approximation of r2.  *)
let rec result_matches r1 r2 =
  match r1, r2 with
  | Any, _ -> true
  | _, Any -> raise (Invalid_argument "result_matches: wild card on rhs")
  | Type t, _ -> result_type_matches t r2
  | Unit, Unit -> true
  | Int i1, Int i2 -> i1 = i2
  | Float f1, Float f2 -> abs_float (f1 -. f2) < !cmp_eps
  | Aggr l1, Aggr l2 -> for_all2_and_same_len result_matches l1 l2
  | Carray v1, Carray v2 -> carrays_equal v1 v2
  | SolverResult r1, SolverResult r2 -> r1 = r2
  | Exn e1, Exn e2 -> exns_equal e1 e2
  | RootInfo r1, RootInfo r2 -> r1 = r2
  | _, _ -> false
and result_type_matches r1 r2 =
  match r1, r2 with
  | Any, _ -> true
  | Type t, _ -> raise (Invalid_argument "result_matches: nested Type")
  | Unit, Unit -> true
  | Int _, Int _ -> true
  | Float _, Float _ -> true
  | Aggr l1, Aggr l2 -> for_all2_and_same_len result_type_matches l1 l2
  | Carray _, Carray _ -> true
  | SolverResult _, SolverResult _ -> true
  | RootInfo _, RootInfo _ -> true
  | Exn e1, Exn e2 -> true
  | _, Any | _, Type _ ->
    raise (Invalid_argument "result_matches: wild card on rhs")
  | _, _ -> false
and for_all2_and_same_len f r1 r2 =
  match r1, r2 with
  | [], [] -> true
  | x::xs, y::ys -> f x y && for_all2_and_same_len f xs ys
  | _, _ -> false
and carrays_equal v1 v2 =
  let n = Carray.length v1 in
  let rec go i =
    if i < n then abs_float (v1.{i} -. v2.{i}) < !cmp_eps && go (i+1)
    else true
  in
  n = Carray.length v2 && go 0
and exns_equal =
  (* NB: exn objects passed through Marshal does not seem to match a pattern
     that it should otherwise match.  For example, an Invalid_argument "foo"
     that came through Marshal does not match the pattern Invalid_argument _.
     This is partly why the generated, compiled code does everything from
     executing the script to comparing results.  *)

  (* Compare only the tags, except for known, common messages.  *)
  (* The table construction is delayed to avoid having this initialization run
     every time the test code starts up.  *)
  let fixed_msgs =
    lazy (let fixed_msgs = Hashtbl.create 10 in
          Hashtbl.add fixed_msgs "index out of bounds" ();
          Hashtbl.add fixed_msgs "hd" ();
          Hashtbl.add fixed_msgs "tl" ();
          Hashtbl.add fixed_msgs "Array.make" ();
          Hashtbl.add fixed_msgs "Bigarray.create: negative dimension" ();
          fixed_msgs)
  in
  fun e1 e2 ->
  match e1, e2 with
  | Failure m1, Failure m2
  | Invalid_argument m1, Invalid_argument m2 ->
    if Hashtbl.mem (Lazy.force fixed_msgs) m1 then m1 = m2
    else if Hashtbl.mem (Lazy.force fixed_msgs) m2 then false
    else true
  | _, _ -> e1 = e2
and cmp_eps = ref 1e-5

(* The quickcheck_main function below is called whenever any of the modules
   Quickcheck_ida_serial, Quickcheck_cvode_serial, etc. is loaded, which serves
   as the program's entry point.  The problem is, this happens even when
   loading from the top level.  This is not only very time consuming, but
   quickcheck_main almost always fails because the top level doesn't set up
   command line arguments properly.  So the loading of the module fails
   altogether.

   To avoid this problem, the following flag inhibits quickcheck_main, turning
   it into a no-op.  The uninhibited quickcheck_main can still be called
   explicitly as quickcheck_main_internal if needed.  *)
let inhibit_quickcheck_main = ref false

(* Model of an imperative language that involves only straight-line code.  Each
   sundials module (CVODE, IDA, etc) should implement this module type and pass
   it to TestImperativeLang below.  *)
module type ImperativeLang =
sig
  type model                (* model of the state of execution environment *)
  type cmd                  (* a command *)

  type script = model * cmd list

  val pp_model : model pp
  val pp_cmd : cmd pp
  val copy_model : model -> model

  val gen_script : unit -> script

  val shrink_script : script -> script Fstream.t

  (* Execute a command and mutate the model.  *)
  val exec_cmd : model -> cmd -> result

  (* Summarize the initial model (return e.g. initial values).  *)
  val model_init : model -> result

  (* The type of a function that compiles a script to OCaml and prints that to
     a file, which is then compiled.  Some functions in TestImperativeLang take
     a function of this type as input.  Ideally, this should really be

       val ml_file_of_script : script -> string -> unit

     but an implementation of that signature generally depends on camlp4, which
     is big enough to slow down linking.  TestImperativeLang is used in both
     the test case generator and the generated test cases; as test cases are
     compiled and linked many times, we want to confine dependence on camlp4 to
     the test generator only.  *)
  type ml_file_of_script = script -> string -> unit

  (* Checks if results match.  The first argument is the expected result, and
     may be a pattern covering several possible results rather than a single
     result.  Float values should be compared with a somewhat large
     tolerance.  *)
  val result_matches : result -> result -> bool
end

(* Functions for generating, compiling, executing, and shrinking test cases of
   ImperativeLang, taking care of the relevant process handling.  *)
module TestImperativeLang (L : ImperativeLang) =
struct
  open L

  let verbose = ref false

  (* Print results in the form
     init:   foo
     step 1: bar
     step 2: baz
     ...
     unless read_write_invariance is set.
  *)
  let pp_results, dump_results, show_results, display_results,
      print_results, prerr_results =
    printers_of_pp (fun ?(prec=0) fmt rs ->
        if !read_write_invariance then pp_vlist pp_result fmt rs
        else
          (* List one result per line, with step numbers starting from 0; as an
             exception, the 0th step is printed "init" rather than "step 0". *)
          let nsteps = List.length rs - 1 in
          let step_width = String.length (string_of_int (nsteps-1)) in
          let pad_show n = let s = string_of_int n in
            String.make (step_width - String.length s) ' ' ^ s
          in
          pp_seq "@[<v>" ";" "@]" fmt
            (Fstream.mapi (fun i r fmt ->
                 if i = 0 then
                   Format.fprintf fmt "init:%s  " (String.make step_width ' ')
                 else
                   Format.fprintf fmt "step %s: " (pad_show i);
                 pp_result fmt r)
                (Fstream.of_list rs)))

  (* Print commands in the form
     step 1: bar
     step 2: baz
     ...
     unless read_write_invariance is set.
  *)
  let pp_cmds, dump_cmds, show_cmds, display_cmds, print_cmds, prerr_cmds =
    printers_of_pp (fun ?(prec=0) fmt cmds ->
        if !read_write_invariance then pp_list pp_cmd fmt cmds
        else
          (* List one command per line, with step numbers starting from 1.  *)
          let nsteps = List.length cmds in
          let step_width = String.length (string_of_int nsteps) in
          let pad_show n = let s = string_of_int n in
            String.make (step_width - String.length s) ' ' ^ s
          in
          pp_seq "@[<v>" ";" "@]" fmt
            (Fstream.mapi (fun i cmd fmt ->
                 Format.fprintf fmt "step %s: " (pad_show (i+1));
                 pp_cmd fmt cmd)
                (Fstream.of_list cmds)))


  (* Print scripts (= model * cmd list) in the form
     model = { ... }
     step 1: foo (* first cmd *)
     step 2: bar (* second cmd *)
     ...
  *)
  let pp_script, dump_script, show_script, display_script,
      print_script, prerr_script
    =
    printers_of_pp (fun ?prec fmt (model, cmds) ->
        Format.fprintf fmt "@[<hov 2>%smodel =@ "
          (if !read_write_invariance
           then "let "
           else "");
        pp_model fmt model;
        Format.fprintf fmt "@]@\n@[<hov 2>%scmds =@ "
          (if !read_write_invariance
           then "let "
           else "");
        pp_cmds fmt cmds;
        if !read_write_invariance
        then Format.fprintf fmt "@]@\n@[<hov 2>let script =@ (model, cmds)@]"
        else Format.fprintf fmt "@]"
      )

  (* Run a list of commands on the model.  *)
  let model_run (model, cmds) =
    let model = copy_model model in
    (* Evaluation order requires this let here.  *)
    let head = model_init model in
    Fstream.cons head (Fstream.map (exec_cmd model) (Fstream.of_list cmds))

  let copy_file ~from_file ~to_file () =
    let infile  = open_in from_file
    and outfile = open_out to_file in
    (try while true do output_char outfile (input_char infile) done
     with End_of_file -> ());
    close_in infile;
    close_out outfile

  (* The main UI to be called from each generated test case.  *)
  let test_case_driver model cmds =
    let just_cmp = ref false
    and orig_model = copy_model model
    and buffer_errors = ref false
    and error_buffer = Buffer.create 1
    in
    Arg.parse
      [("--just-cmp", Arg.Set just_cmp,
        "don't print anything; just give the exit code");
       ("--buffer-errors", Arg.Set buffer_errors,
        "for internal use only");
       ("--read-write-invariance", Arg.Set read_write_invariance,
        "print data in a format that can be fed to ocaml toplevel")]
      (fun _ -> ()) "a test case generated by quickcheck";
    let step = ref 0 in
    let mismatches = ref [] in
    let expected_results = ref [] in
    let actual_results = ref [] in
    let errout =
      if !buffer_errors then Format.formatter_of_buffer error_buffer
      else Format.err_formatter
    in
    let do_cmd thunk =
      let expected = if !step = 0 then model_init model
        else exec_cmd model cmds.(!step - 1)
      and actual = try Lazy.force thunk with exn -> Exn exn in
      if !just_cmp then
        (if not (result_matches expected actual)
         then exit 1)
      else
        (actual_results := actual::!actual_results;
         expected_results := expected::!expected_results;
         if not (result_matches expected actual)
         then mismatches := (!step, expected, actual)::!mismatches);
      step := !step + 1
    in
    let finish () =
      if not !just_cmp then
        match !mismatches with
        | [] -> Printf.printf "OK, test successful.\n"; 0
        | ms ->
          let prerrf fmt = Format.fprintf Format.err_formatter fmt in
          let prerr_mismatch (i, exp, act) =
            prerrf "Result mismatch on step %d:@\ngot@\n  %s@\n\
                    but expected@\n  %s@\n"
              i (show_result act) (show_result exp)
          in
          Format.pp_print_flush errout ();
          prerrf "Test failed.@\n[Reason]@\n";
          List.iter prerr_mismatch !mismatches;
          prerrf "@\n[Test Case]@\n";
          prerr_script (orig_model, Array.to_list cmds);
          prerrf "@\n@\n[Program Output]@\n";
          (* If we print errors straight to standard error, they will appear
             before the "[Test Case]" above, which makes it easy to miss.
             If at all possible, we buffer the output and put it here, under
             the [Program Output] section.  *)
          if !buffer_errors then
            prerrf "%s" (Buffer.contents error_buffer);
          prerr_results (List.rev !actual_results);
          prerrf "@\n@\n[Expected Output]@\n";
          prerr_results (List.rev !expected_results);
          prerrf "@\n";
          1
      else 0
    in
    let err_handler details =
      (* Follows the format of sundials 2.5.0.  *)
      Format.fprintf errout "\n[%s ERROR]  %s\n  %s\n\n"
        details.Ida.module_name details.Ida.function_name
        details.Ida.error_message
    in
    do_cmd, finish, err_handler

    (* Close file descriptor when f finishes, whether normally or by an
       exception.  *)
    let with_file_descr descr f =
      let ret = try f descr with exn -> Unix.close descr; raise exn in
      Unix.close descr;
      ret

    (* Code generation options *)
    let test_exec_file = ref "./tmp"
    let test_compiler = ref "false"
    let test_failed_file = ref "./failed.ml"

    let compile mk_ml_file =
      let test_src_file = !test_exec_file ^ ".ml" in
      mk_ml_file test_src_file;
      let compile_command =
        !test_compiler ^ " -o " ^ !test_exec_file ^ " " ^ test_src_file
      in
      match Unix.system compile_command with
      | Unix.WEXITED 0 -> ()
      | Unix.WEXITED _ | Unix.WSIGNALED _ | Unix.WSTOPPED _ ->
        copy_file ~from_file:test_src_file ~to_file:!test_failed_file ();
        raise (Failure ("Internal error: generated code failed to compile.  \
                         Compilation command was: " ^ compile_command))

    (* Just checks if the generated code's self-test is successful.  *)
    let prop_script_ok (ml_file_of_script : ml_file_of_script) script =
      (try compile (ml_file_of_script script)
       with exn -> raise (AbortTests exn));
      let dev_null =
        try Unix.openfile "/dev/null" [Unix.O_RDWR; Unix.O_TRUNC] 0
        with Unix.Unix_error (Unix.ENOENT, "open", "/dev/null") ->
          (* Windows *)
          try Unix.openfile "NUL" [Unix.O_RDWR; Unix.O_TRUNC] 0
          with Unix.Unix_error (Unix.ENOENT, "open", "NUL") ->
            Printf.fprintf stderr
              "ERROR: Can't find /dev/null or equivalent on your system.  \
               Giving up.  Fix prop_ida_ok in quickcheck_ida.ml and try again.";
            exit 2
      in
      with_file_descr dev_null (fun dev_null ->
          let pid =
            Unix.create_process !test_exec_file
              [|!test_exec_file; "--just-cmp"|]
              dev_null dev_null dev_null
          in
          let _, exit_code = Unix.waitpid [] pid in
          if exit_code = Unix.WEXITED 0 then OK
          else Falsified exit_code)

    let quickcheck_scripts ml_file_of_script shrink max_tests =
      let err = Format.err_formatter in
      let fprintf = Format.fprintf err in
      let prop = prop_script_ok ml_file_of_script in
      let shrinker = if shrink then shrink_script else (fun _ -> Fstream.nil)
      in
      let result =
        if !verbose
        then quickcheck gen_script shrinker ~pp_input:pp_script prop max_tests
        else quickcheck gen_script shrinker prop max_tests
      in
      match result with
      | None | Some (_, OK) -> None
      | Some (_, Failed exn) ->
        fprintf "internal error; uncaught exception in property";
        raise exn
      | Some (script, Falsified (Unix.WEXITED 1)) ->
        (* Exit code 1 means the test ran to completion but there was a mismatching
           result.  We can let the test code print the report.  Note we need to
           re-generate the code here because the failed script is generally not the
           last one that was tried, and its source code may therefore have been
           overwritten.  *)
        compile (ml_file_of_script script);
        copy_file ~from_file:(!test_exec_file ^ ".ml")
          ~to_file:!test_failed_file ();
        let exit_status =
          Unix.system
            (Printf.sprintf "%s --buffer-errors %s"
               !test_exec_file
               (if !read_write_invariance then " --read-write-invariance"
                else ""))
        in
        assert (exit_status = Unix.WEXITED 1);
        Format.fprintf err
          "@\n@\nFailed test code saved in %s, compile with:@\n%s@\n"
          !test_failed_file (!test_compiler ^ " " ^ !test_failed_file);
        Some script
      | Some (script, Falsified stat) ->
        (* Other exit statuses mean the test case crashed.  We'll run it again,
           this time letting it produce as much output as it can, and compare that
           to the expected results.  *)
        compile (ml_file_of_script script);
        copy_file ~from_file:(!test_exec_file ^ ".ml")
          ~to_file:!test_failed_file
          ();
        Format.fprintf err "@\n[Reason]@\nTest code %s.@\n\
                            @\n[Test Case]@\n"
          (match stat with
           | Unix.WEXITED n -> Printf.sprintf "crashed (exit code %d)" n
           | Unix.WSIGNALED n -> Printf.sprintf "received signal %d" n
           | Unix.WSTOPPED n -> Printf.sprintf "stopped by signal %d" n);
        pp_script err script;
        Format.fprintf err "@\n@\n[Program Output]@\n";
        Format.pp_print_flush err ();
        flush stderr;
        (* No --buffer-errors here, because the test code will likely not get to
           the point where it flushes the buffer.  *)
        ignore (Unix.system
                  (if !read_write_invariance
                   then !test_exec_file ^ " --read-write-invariance"
                   else !test_exec_file));
        flush stdout;
        flush stderr;
        Format.fprintf err "@\n[Expected Output]@\n";
        pp_results err (Fstream.to_list (model_run script));
        Format.pp_print_newline err ();
        Format.fprintf err
          "@\n@\nFailed test code saved in %s, compile with:@\n%s@\n"
          !test_failed_file (!test_compiler ^ " " ^ !test_failed_file);
        Some script


    (* Entry point for the test generator.  ml_file_of_script receives the
       random seed as the first argument, for informative purposes.  *)
    let quickcheck_main_internal ml_file_of_script () =
      let randseed =
        Random.self_init ();
        ref (Random.int ((1 lsl 30) - 1))
      and max_tests = ref 50
      and shrink = ref true
      in
      let options = [("--exec-file", Arg.Set_string test_exec_file,
                      "test executable name \
                       (must be absolute, prefixed with ./, or on path)");
                     ("--no-shrink", Arg.Clear shrink,
                      "don't shrink the failed test case");
                     ("--failed-file", Arg.Set_string test_failed_file,
                      "file in which to dump the failed test case");
                     ("--compiler", Arg.Set_string test_compiler,
                      "compiler name with compilation options");
                     ("--rand-seed", Arg.Set_int randseed,
                      "seed value for random generator");
                     ("--verbose", Arg.Set verbose,
                      "print each test script before trying it");
                     ("--read-write-invariance", Arg.Set read_write_invariance,
                      "print data in a format that can be fed to ocaml toplevel");
                    ] in
      Arg.parse options (fun n -> max_tests := int_of_string n)
        "randomly generate programs using CVODE and check if they work as expected";

      Printf.printf "random generator seed value = %d\n" !randseed;
      flush stdout;
      Random.init !randseed;
      size := 1;

      let ml_file_of_script script file =
        ml_file_of_script script file;
        let chan = open_out_gen [Open_text; Open_append; Open_wronly] 0 file in
        Printf.fprintf chan
          "\n(* generated with random seed %d, test case %d *)\n"
          !randseed !test_case_number;
        close_out chan
      in

      (* If loaded from top level,  *)
      quickcheck_scripts ml_file_of_script !shrink !max_tests

    let quickcheck_main ml_file_of_script () =
      if !inhibit_quickcheck_main then None
      else quickcheck_main_internal ml_file_of_script ()
  end
