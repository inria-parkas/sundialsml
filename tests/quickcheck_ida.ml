module Carray = Sundials.Carray
open Pprint
open Quickcheck
open Quickcheck_sundials
module Roots = Sundials.Roots
module RootDirs = Sundials.RootDirs

(* A state-machine model for IDA sessions.  *)
type session_model =
  {
    resfn : resfn_type;
    mutable solver : Ida.linear_solver;
    mutable solving : bool;
    mutable consistent : bool;
    mutable last_query_time : float;
    mutable last_tret : float;
    mutable roots : (float * Ida.Roots.root_event) array;
    mutable root_info : Roots.t;
    vec : Carray.t;
    vec' : Carray.t;
    t0 : float;
    vec0 : Carray.t;
    vec'0 : Carray.t;
  }
and script = session_model * cmd list
and cmd = SolveNormal of float
          | GetRootInfo
and result = Unit | Int of int | Float of float
              | Any
              | Type of result
              | Aggr of result list
              | Carray of Carray.t    (* NB: always copy the array! *)
              | SolverResult of Ida.solver_result
              | Exn of exn
              | RootInfo of Ida.Roots.t
and resfn_type = ResFnLinear of Carray.t

let carray x = Carray (Carray.of_carray x)

  (* Whole-test results *)
type failure_type = ResultMismatch of int * result * result
                    | TestCodeDied of result list
                    | TestCodeOverrun

(* Model *)

let gen_resfn neqs () =
  gen_choice
    [|
      fun () -> let v = Carray.create neqs in
                for i = 0 to Carray.length v - 1 do
                  v.{i} <- float_of_int i
                done;
                ResFnLinear v
    |]
    ()

let gen_solver lapack neqs =
  match Random.int 1 with
  | 0 -> if lapack && Random.bool () then Ida.LapackDense else Ida.Dense
  | _ -> assert false

let gen_roots () =
  (* We have to ensure the roots don't coincide with any of the query times or
     with one another, for otherwise the order in which they fire is dictated
     by floating point error and is unpredictable.  *)
  let gen_root () =
    (gen_root_time (),
     gen_choice [| Roots.Rising; Roots.Falling |])
  in
  match Random.int 3 with
  | 0 -> [||]
  | _ -> uniq_array (gen_array gen_root ())

let init_vec_for neqs t0 resfn =
  match resfn with
  | ResFnLinear slopes -> Carray.init neqs 0., Carray.of_carray slopes

let gen_model () =
  let neqs  = min 10 (gen_pos ()) in
  let resfn = gen_resfn neqs () in
  let t0 = gen_t0 () in
  let vec0, vec'0  = init_vec_for neqs t0 resfn in
  let neqs  = Carray.length vec0 in
  let roots = gen_roots () in
  {
    resfn = resfn;
    solver = gen_solver false neqs;
    solving = false;
    last_query_time = t0;
    last_tret = t0;
    consistent = true;
    roots = roots;
    root_info = Roots.create (Array.length roots);
    vec   = Carray.of_carray vec0;
    vec'  = Carray.of_carray vec'0;
    t0 = t0;
    vec0  = vec0;
    vec'0 = vec'0;
  }

let shrink_root_dir = function
  | Roots.Rising -> Fstream.nil
  | Roots.Falling -> Fstream.singleton Roots.Rising
  | Roots.NoRoot -> assert false
let shrink_root_pos = shrink_root_time

let shrink_model model cmds =
  let update_roots model roots =
    let root_info = Roots.create (Array.length roots) in
    { model with roots = roots; root_info = root_info }
  in
  Fstream.map (fun t -> ({ model with last_tret = t; t0 = t; }, cmds))
    (shrink_discrete_float model.t0)
  @@ Fstream.map (fun roots -> (update_roots model roots, cmds))
    (shrink_array (shrink_pair shrink_root_pos shrink_root_dir)
       model.roots)

let copy_resfn = function
  | ResFnLinear slopes -> ResFnLinear (Carray.of_carray slopes)
let copy_model m =
  {
    resfn = copy_resfn m.resfn;
    solver = m.solver;
    solving = m.solving;
    last_query_time = m.last_query_time;
    last_tret = m.last_tret;
    roots = Array.copy m.roots;
    root_info = Roots.copy m.root_info;
    consistent = m.consistent;
    vec = Carray.of_carray m.vec;
    vec' = Carray.of_carray m.vec';
    t0 = m.t0;
    vec0 = Carray.of_carray m.vec0;
    vec'0 = Carray.of_carray m.vec'0;
  }

(* Commands *)

let gen_cmd =
  let cases =
    [| (fun () -> SolveNormal (gen_query_time ()));
       (fun () -> GetRootInfo) |]
  in
  fun () -> gen_choice cases ()

let shrink_cmd = function
  | SolveNormal dt ->
    Fstream.map (fun dt -> SolveNormal dt) (shrink_query_time dt)
  | GetRootInfo -> Fstream.nil

let gen_cmds = gen_list gen_cmd

(* Scripts (model + command list) *)

let gen_script () =
  let model = gen_model ()
  and cmds  = gen_cmds ()
  in (model, cmds)

let shrink_neqs model cmds =
  (* Reduce commands that are dependent on number of equations.  *)
  let drop_from_cmd i = function
    | cmd -> cmd
  in
  let copy_vec_drop i v =
    let n = Carray.length v in
    let w = Carray.create (n-1) in
    for j = 0 to i-1 do
      w.{j} <- v.{j}
    done;
    for j = i+1 to n - 1 do
      w.{j-1} <- v.{j}
    done;
    w
  in
  let copy_resfn_drop i = function
    | ResFnLinear slopes -> ResFnLinear (copy_vec_drop i slopes)
  in
  let drop_eq i =
    let model =
      {
        resfn = copy_resfn_drop i model.resfn;
        roots = Array.copy model.roots;
        root_info = Roots.copy model.root_info;
        solving = model.solving;
        consistent = model.consistent;
        vec = copy_vec_drop i model.vec;
        vec' = copy_vec_drop i model.vec';
        t0 = model.t0;
        vec0 = copy_vec_drop i model.vec0;
        vec'0 = copy_vec_drop i model.vec'0;
        solver = model.solver;
        last_query_time = model.last_query_time;
        last_tret = model.last_tret;
      }
    in (model, List.map (drop_from_cmd i) cmds)
  in
  let n = Carray.length model.vec in
  Fstream.guard (n > 1) (Fstream.map drop_eq (Fstream.enum 0 (n - 1)))

let shrink_script (model, cmds) =
  shrink_neqs model cmds
  @@ shrink_model model cmds
  @@ Fstream.map (fun cmds -> (model, cmds)) (shrink_list shrink_cmd cmds)

(* Pretty-printing and result comparison.  *)

let cmp_eps = ref 1e-5

let pp_solver, dump_solver, show_solver, display_solver,
  print_solver, prerr_solver =
  printers_of_show (fun s ->
    let solver_name = function
      | Ida.Dense -> "Dense"
      | Ida.Band range -> Printf.sprintf "Band { mupper=%d; mlower=%d }"
        range.Ida.mupper range.Ida.mlower
      | Ida.Sptfqmr _ | Ida.Spbcg _ | Ida.Spgmr _
      | Ida.LapackBand _ | Ida.LapackDense _ ->
        raise (Failure "linear solver not implemented")
    in
    if !read_write_invariance then "Ida." ^ solver_name s
    else solver_name s)

let pp_ida_ident fmt ident =
  if !read_write_invariance then Format.fprintf fmt "Ida.%s" ident
  else Format.fprintf fmt "%s" ident

let pp_result, dump_result, show_result, display_result,
  print_result, prerr_result =
  let rec pre_pp_result arg_pos fmt = function
    | Any -> Format.fprintf fmt "_"
    | Unit -> Format.pp_print_string fmt (if !read_write_invariance
                                          then "Unit"
                                          else "()")
    | Int i -> pp_parens (arg_pos && i < 0) pp_int fmt i
    | Float f -> pp_parens (arg_pos && f < 0.) pp_float fmt f
    | Type r -> pp_parens arg_pos (fun fmt r ->
      pp_ida_ident fmt "Type ";
      pre_pp_result true fmt r) fmt r
    | Carray ca -> pp_carray fmt ca
    | SolverResult Ida.Continue -> pp_ida_ident fmt "Continue"
    | SolverResult Ida.RootsFound -> pp_ida_ident fmt "RootsFound"
    | SolverResult Ida.StopTimeReached -> pp_ida_ident fmt "StopTimeReached"
    | RootInfo roots -> pp_parens arg_pos (fun fmt roots ->
      pp_string_verbatim fmt "RootInfo ";
      pp_root_info fmt roots) fmt roots
    | Aggr rs -> pp_parens arg_pos (fun fmt rs ->
      pp_string_verbatim fmt "Aggr ";
      pp_list (pre_pp_result false) fmt rs) fmt rs
    | Exn exn -> pp_parens arg_pos (fun fmt exn ->
      pp_string_verbatim fmt "exception ";
      pp_string_verbatim fmt (Printexc.to_string exn)) fmt exn
  in printers_of_pp (pre_pp_result false)
let pp_results, dump_results, show_results, display_results,
  print_results, prerr_results =
  printers_of_pp (pp_list pp_result)

let pp_resfn_type, dump_resfn_type, show_resfn_type, display_resfn_type,
  print_resfn_type, prerr_resfn_type
    =
  printers_of_pp
    (fun fmt -> function
    | ResFnLinear slope -> pp_string_verbatim fmt "ResFnLinear ";
      pp_carray fmt slope
    )

let pp_cmd, dump_cmd, show_cmd, display_cmd, print_cmd, prerr_cmd =
  printers_of_pp (fun fmt -> function
  | SolveNormal f ->
    Format.fprintf fmt "SolveNormal ";
    if f < 0. then Format.fprintf fmt "(";
    pp_float fmt f;
    if f < 0. then Format.fprintf fmt ")"
  | GetRootInfo -> Format.fprintf fmt "GetRootInfo")

let pp_cmds, dump_cmds, show_cmds, display_cmds, print_cmds, prerr_cmds =
  printers_of_pp (fun fmt cmds ->
    if !read_write_invariance then pp_list pp_cmd fmt cmds
    else
      (* List one command per line, with step numbers starting from 0.  *)
      let nsteps = List.length cmds in
      let step_width = String.length (string_of_int (nsteps - 1)) in
      let pad_show n = let s = string_of_int n in
                       String.make (step_width - String.length s) ' ' ^ s
      in
      pp_seq "[" ";" "]" fmt
        (Fstream.mapi (fun i cmd fmt ->
          Format.fprintf fmt "Step %s: " (pad_show i);
          pp_cmd fmt cmd)
           (Fstream.of_list cmds)))

let pp_model, dump_model, show_model, display_model,
    print_model, prerr_model
      =
  printers_of_pp
    (let descriptive_fields =
       ["resfn", (fun fmt m -> pp_of_show show_resfn_type fmt m.resfn);
        "solver", (fun fmt m -> pp_of_show show_solver fmt m.solver);
        "roots", (fun fmt m -> pp_array (pp_pair pp_float pp_root_event)
                                 fmt m.roots);
        "vec", (fun fmt m -> pp_carray fmt m.vec);
        "vec'", (fun fmt m -> pp_carray fmt m.vec');
        "t0", (fun fmt m -> pp_float fmt m.t0);
        "vec0", (fun fmt m -> pp_carray fmt m.vec0);
        "vec'0", (fun fmt m -> pp_carray fmt m.vec'0);
       ]
     and dynamic_fields =
       ["root_info", (fun fmt m -> pp_root_info fmt m.root_info);
        "last_query_time", (fun fmt m -> pp_float fmt m.last_query_time);
        "last_tret", (fun fmt m -> pp_float fmt m.last_tret);
        "solving", (fun fmt m -> pp_bool fmt m.solving);
        "consistent", (fun fmt m -> pp_bool fmt m.consistent);
       ]
     in
     let rw_invar = pp_record (descriptive_fields @ dynamic_fields)
     and concise = pp_record descriptive_fields in
     fun fmt m ->
       if !read_write_invariance then rw_invar fmt m
       else concise fmt m)

let pp_script, dump_script, show_script, display_script,
    print_script, prerr_script
      =
  printers_of_pp (fun fmt (model, cmds) ->
    Format.fprintf fmt "@[<hov 2>let model =@ ";
    pp_model fmt model;
    Format.fprintf fmt "@]@\n@[<hov 2>let cmds =@ ";
    pp_cmds fmt cmds;
    Format.fprintf fmt "@]@\n@[<hov 2>let script =@ (model, cmds)@]"
  )

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
  | Exn e1, Exn e2 -> raise (Invalid_argument "result_matches: Type Exn")
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
and exns_equal e1 e2 =
  (* Compare only the tags *)
  match e1, e2 with
  | Failure _, Failure _ -> true
  | Invalid_argument _, Invalid_argument _ -> true
  | _, _ -> e1 = e2

let is_exn = function
  | Exn _ -> true
  | _ -> false
let not_exn x = not (is_exn x)

(* Compilation and execution *)

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
    raise (Failure ("Internal error: generated code failed to compile.  \
                     Compilation command was: " ^ compile_command))

(* Compile test code, run it, and return a list of responses.  *)
let compile_run ml_code =
  compile ml_code;
  let (pipe_stdout, pipe_stdin, pipe_stderr) as pipes =
    Unix.open_process_full
      (!test_exec_file ^ " --marshal-results")
      (Unix.environment ())
  in
  let recv () =
    try Some (Marshal.from_channel pipe_stdout : result)
    with End_of_file -> (ignore (Unix.close_process_full pipes); None)
  in
  Fstream.generate recv

(* Model interpretation *)

let exact_soln typ t0 vec0 vec'0 t vec vec' =
  match typ with
  | ResFnLinear slopes ->
    for i = 0 to Carray.length vec - 1 do
      vec'.{i} <- slopes.{i};
      vec.{i} <- slopes.{i} *. (t -. t0);
    done
let exact_init model =
  match model.resfn with
  | ResFnLinear slopes ->
    Carray.fill model.vec 0.;
    Carray.fill model.vec0 0.;
    model.consistent <- true

(* Find the smallest root(s) in the range (t1..t2].  *)
let find_roots roots t1 t2 =
  let n = Array.length roots in
  let pos i = fst roots.(i) in
  if n = 0 then []
  else
    let record  = ref t2
    and holders = ref []
    in
    for i = 0 to n-1 do
      if t1 < pos i && pos i <= !record
      then
        begin
          if pos i < !record then
            (record := pos i;
             holders := [i])
          else
            holders := i::!holders
        end
    done;
    !holders

(* Run a single command on the model.  *)
let model_cmd model = function
  | SolveNormal dt ->
    (* NB: we don't model interpolation failures -- t will be monotonically
       increasing.  *)
    let t = model.last_query_time +. dt in
    (* Undocumented behavior (sundials 2.5.0): solve_normal with t=t0 usually
       fails with "tout too close to t0 to start integration", but sometimes
       succeeds.  Whether it succeeds seems to be unpredictable.  *)
    if t = model.t0 then Any
    else
      let tret, flag =
        (* NB: roots that the solver is already sitting on are ignored, unless
           dt = 0.  *)
        let tstart =
          if dt = 0.
          then model.last_tret -. time_epsilon
          else model.last_tret
        in
        match find_roots model.roots tstart t with
        | [] ->
          (* Undocumented behavior (sundials 2.5.0): a non-root return from
             solve_normal resets root info to Rising.  *)
          Roots.fill_all model.root_info Roots.Rising;
          t, SolverResult Ida.Continue
        | (i::_) as is ->
          let tret = fst model.roots.(i) in
          Roots.reset model.root_info;
          List.iter
            (fun i -> Roots.set model.root_info i (snd model.roots.(i)))
            is;
          (tret,
           (* If the queried time coincides with a root, then it's reasonable
              for the solver to prioritize either.  In real code, this is most
              likely determined by the state of the floating point error.  *)
           if tret = t then Type (SolverResult Ida.Continue)
           else SolverResult Ida.RootsFound)
      in
      model.last_query_time <- t;
      model.last_tret <- tret;
      model.solving <- true;
      exact_soln model.resfn model.t0 model.vec0 model.vec'0
                                 tret model.vec  model.vec';
      Aggr [Float tret; flag; carray model.vec; carray model.vec']
  | GetRootInfo ->
    if model.solving then RootInfo (Roots.copy model.root_info)
    else
      (* FIXME: this should be Exn Ida.IllInput or something like that.  *)
      Type (RootInfo (Roots.copy model.root_info))

(* Run a list of commands on the model.  *)
let model_run (model, cmds) =
  let model = copy_model model in
  (* Evaluation order requires this let here.  *)
  let head = Aggr [Float model.t0;
                   carray model.vec;
                   carray model.vec']
  in
  Fstream.cons head (Fstream.map (model_cmd model) (Fstream.of_list cmds))

(* Checking results *)
let compare_results model_run compiled_run =
  let rec cmp i ms cs =
    match Fstream.decons ms, Fstream.decons cs with
    | None, None -> OK
    | Some (m, ms), Some (c, cs) ->
      if result_matches m c then cmp (i+1) ms cs
      else Falsified (ResultMismatch (i, m, c))
    | Some _, None ->
      Falsified (TestCodeDied (Fstream.to_list ms))
    | None, Some _ ->
      Falsified TestCodeOverrun
  in
  cmp 0 model_run compiled_run

let verbose = ref false

let quickcheck_ida ml_file_of_script max_tests =
  let prop script =
    compare_results (model_run script) (compile_run (ml_file_of_script script))
  in
  let err = Format.err_formatter in
  let fprintf = Format.fprintf err in
  let result =
    if !verbose then
      quickcheck gen_script shrink_script ~pp_input:pp_script prop max_tests
    else
      quickcheck gen_script shrink_script prop max_tests
  in
  match result with
  | None | Some (_, OK) -> result
  | Some (_, Failed exn) ->
    fprintf "internal error; uncaught exception in property";
    raise exn
  | Some (script, Falsified reason) ->
    Format.fprintf err "Failed test code saved in %s.\n\n[Reason]\n"
      !test_failed_file;
    (match reason with
    | TestCodeDied rs ->
      Format.fprintf err "Test code exited unexpectedly.  It was \
                          supposed to produce the following additional \
                          results:\n";
      pp_results err rs
    | TestCodeOverrun -> assert false   (* Shouldn't happen.  *)
    | ResultMismatch (i,e,a) ->
      Format.fprintf err "Result mismatch on %s:\nexpected\n%s\nbut got\n%s\n"
        (if i = 0 then "init" else ("step " ^ string_of_int i))
        (show_result e)
        (show_result a));
    Format.pp_print_string err "\n[Test Case]\n";
    pp_script err script;
    Format.pp_print_string err "\n\n[Program Output]\n";
    Format.pp_print_flush err ();
    (* The test file contains the last script tried, not the last script that
       failed.  We need to reinstate the failing script before running it again
       to retrieve its output.  *)
    compile (ml_file_of_script script);
    let pid = Unix.create_process !test_exec_file
      (if !read_write_invariance
       then [|!test_exec_file; "--read-write-invariance"|]
       else [|!test_exec_file|])
      Unix.stdin Unix.stdout Unix.stderr
    in
    (* Copy the failed test from test_exec_file.ml to test_failed_file.  *)
    let _ = Unix.waitpid [] pid in
    let infile  = open_in (!test_exec_file ^ ".ml")
    and outfile = open_out !test_failed_file in
    let _ =
      try while true do output_char outfile (input_char infile) done
      with End_of_file -> ()
    in
    flush stderr;
    Format.pp_print_string err "\n\n[Expected Output]\n";
    pp_results err (Fstream.to_list (model_run script));
    Format.pp_print_flush err ();
    result
