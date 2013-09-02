module Carray = Sundials.Carray
open Pprint
open Quickcheck
open Quickcheck_sundials
module Roots = Sundials.Roots
module RootDirs = Sundials.RootDirs

(* Derive pretty-printers for types imported from the Ida module.  It would be
   ideal if we could have the camlp4 extension extract type definitions from
   ../ida.ml and do everything automatically, but that's too much work.  For
   now, we'll make do by manually copying the definitions; the types aren't
   supposed to change often, and the compiler will detect bit rot in this
   case, so it should be OK.  *)
module PPIda = struct
  open Ida (* deriving needs access to type names like linear_solver *)

  (* Copied from ida.mli; must be kept in sync with that file's definition.  *)
  type linear_solver =
    | Dense
    | Band of bandrange
    | LapackDense
    | LapackBand of bandrange
    | Spgmr of sprange
    | Spbcg of sprange
    | Sptfqmr of sprange
  and bandrange = { mupper : int; mlower : int; }
  and sprange = int
  and solver_result =
    | Continue
    | RootsFound
    | StopTimeReached
  and error_details = {
    error_code : int;
    module_name : string;
    function_name : string;
    error_message : string;
  }
  and integrator_stats = {
    num_steps : int;
    num_res_evals : int;
    num_lin_solv_setups : int;
    num_err_test_fails : int;
    last_order : int;
    current_order : int;
    actual_init_step : float;
    last_step : float;
    current_step : float;
    current_time : float
  }
  deriving (pretty ~erase_typedefs:true ~prefix:Ida)
end
open PPIda        (* Bring the printers -- and only the printers into scope. *)

(* A state-machine model for IDA sessions.  *)
type model =
  {
    resfn : resfn_type;
    mutable solver : Ida.linear_solver;

    (* Set when SolveNormal or CalcIC_* has been called on this model.  *)
    mutable solving : bool;

    (* Set when the current state vector satisfies the DAE.  *)
    mutable consistent : bool;
    mutable last_query_time : float;
    mutable last_tret : float;
    mutable roots : (float * Ida.Roots.root_event) array;
    mutable root_dirs : Ida.root_direction array;
    mutable root_info : Roots.t;
    mutable root_info_valid : bool;
    vec : Carray.t;
    vec' : Carray.t;
    t0 : float;
    vec0 : Carray.t;
    vec'0 : Carray.t;

    (* Set during command generation to force the next query time to a
       particular value, if there is a query.  *)
    mutable next_query_time : float option;
  }
and script = model * cmd list
and cmd = SolveNormal of float
        | GetRootInfo
        | SetRootDirection of Ida.root_direction array
        | SetAllRootDirections of Ida.root_direction
        | GetNRoots
        | CalcIC_Y of float * ic_buf
and ic_buf = GetCorrectedIC
           | Don'tGetCorrectedIC
           | GiveBadVector of int (* Pass in a vector of incorrect size. *)
and result = Unit
           | Int of int | Float of float
           | Any
           | Type of result
           | Aggr of result list
           | Carray of Carray.t    (* NB: always copy the array! *)
           | SolverResult of Ida.solver_result
           | Exn of exn
           | RootInfo of Roots.t
and resfn_type =
  | ResFnLinear of Carray.t (* The DAE is
                                 y'.{i} = coefs.{i}
                               with general solution
                                 y.{i} = coefs.{i} * t + C
                               ic_calc_y fails (because y is unused)
                               ic_calc_ya_yd' should work.
                             *)
  | ResFnExpDecay of Carray.t
deriving (pretty ~rename:(Carray.t to carray,
                          Ida.Roots.root_event to root_event,
                          Ida.root_direction to root_direction,
                          Roots.t to root_info,
                          Ida.linear_solver to linear_solver,
                          Ida.solver_result to solver_result)
                 ~optional:(Int, Float, solving, consistent, next_query_time,
                            last_tret, root_info_valid))

let carray x = Carray (Carray.of_carray x)

let copy_resfn = function
  | ResFnLinear slopes -> ResFnLinear (Carray.of_carray slopes)
  | ResFnExpDecay coefs -> ResFnExpDecay (Carray.of_carray coefs)
let copy_model m =
  {
    (* Everything is copied explicitly here instead of using { m with ... }
       because this way the compiler forces us to inspect every new field when
       it's added, so we won't forget to deep copy any of them.  *)
    resfn = copy_resfn m.resfn;
    solver = m.solver;
    solving = m.solving;
    last_query_time = m.last_query_time;
    last_tret = m.last_tret;
    roots = Array.copy m.roots;
    root_dirs = Array.copy m.root_dirs;
    root_info = Roots.copy m.root_info;
    root_info_valid = m.root_info_valid;
    consistent = m.consistent;
    vec = Carray.of_carray m.vec;
    vec' = Carray.of_carray m.vec';
    t0 = m.t0;
    vec0 = Carray.of_carray m.vec0;
    vec'0 = Carray.of_carray m.vec'0;
    next_query_time = m.next_query_time;
  }

let gen_resfn neqs () =
  let _ () =
    (* This is a no-op that causes the compiler to direct you here whenever you
       add a new kind of residual function.  *)
    match ResFnLinear (Carray.create 0) with
    | ResFnLinear _ -> ()
    | ResFnExpDecay _ -> ()
  in
  gen_choice
    [|
      (fun () ->
         let v = Carray.create neqs in
         for i = 0 to Carray.length v - 1 do
           v.{i} <- float_of_int (i+1)
         done;
         ResFnLinear v);
      (fun () ->
         let v = Carray.create neqs in
         (* All coefficients must be non-negative.  *)
         for i = 0 to Carray.length v - 1 do
           v.{i} <- float_of_int (i+1)
         done;
         ResFnExpDecay v);
    |]
    ()

let gen_solver lapack neqs =
  match Random.int 1 with
  | 0 -> if lapack && Random.bool () then Ida.LapackDense else Ida.Dense
  | _ -> assert false

let gen_roots t0 () =
  let gen_root () =
    (gen_root_time t0 (),
     gen_choice [| Roots.Rising; Roots.Falling |])
  in
  match Random.int 3 with
  | 0 -> [||]
  | _ -> uniq_array (gen_array gen_root ())

(* Returns vec0, vec'0 that satisfies the given resfn.  *)
let init_vec_for neqs t0 resfn =
  match resfn with
  | ResFnLinear slopes -> Carray.init neqs 0., Carray.of_carray slopes
  | ResFnExpDecay coefs ->
    let vec0 = Carray.init neqs 1. in
    let vec'0 = Carray.of_carray coefs in
    Carray.mapi (fun i vec0_i -> -. coefs.{i} *. vec0.{i}) vec'0;
    vec0, vec'0

let gen_model () =
  let neqs  = min 10 (gen_pos ()) in
  let resfn = gen_resfn neqs () in
  let t0 = gen_t0 () in
  let vec0, vec'0  = init_vec_for neqs t0 resfn in
  let neqs  = Carray.length vec0 in
  let roots = gen_roots t0 () in
  let num_roots = Array.length roots in
  {
    resfn = resfn;
    solver = gen_solver false neqs;
    solving = false;
    last_query_time = t0;
    last_tret = t0;
    consistent = true;
    roots = roots;
    root_info = Roots.create num_roots;
    root_info_valid = false;
    root_dirs = Array.make num_roots RootDirs.IncreasingOrDecreasing;
    vec   = Carray.of_carray vec0;
    vec'  = Carray.of_carray vec'0;
    t0 = t0;
    vec0  = vec0;
    vec'0 = vec'0;
    next_query_time = None;
  }

(* Commands: note that which commands can be tested without trouble depends on
   the state of the model.  *)

let gen_cmd =
  let cases =
    [| (fun model ->
          let t =
            match model.next_query_time with
            | None -> gen_query_time model.last_query_time ()
            | Some t -> t
          in ({ model with last_query_time = t }, SolveNormal t));
       (fun model -> (model, GetRootInfo));
       (fun model -> (model, GetNRoots));
       (fun model ->
          (* 20% of the time, we (may) generate an incorrectly sized array.  *)
          let size = if Random.int 100 < 20 then gen_nat ()
                     else Array.length model.roots
          in
          let dirs = gen_array ~size:size gen_root_direction ()
          in (model, SetRootDirection dirs));
       (fun model -> (model, SetAllRootDirections (gen_root_direction ())));
       (fun model ->
          match model.resfn with
          (* ic_calc_y doesn't work on ResFnLinear. *)
          | ResFnLinear _ -> (model, GetRootInfo)
          | _ ->
            let t = gen_query_time (model.last_query_time +. discrete_unit) ()
            and rand = Random.int 100 in
            let bad_size () =
              let n = gen_pos () in
              if n <= Carray.length model.vec then n-1
              else n
            in
            let ic_buf =
              (* With 40% probability, don't receive corrected vector.
                 With 10% probability, pass in bad vector.  *)
              if rand < 40 then Don'tGetCorrectedIC
              else if rand < 50 then GiveBadVector (bad_size ())
              else GetCorrectedIC
            in
            (model, CalcIC_Y (t, ic_buf)));
    |]
  in
  fun model -> gen_choice cases model

let shrink_ic_buf model = function
  | GetCorrectedIC -> Fstream.singleton Don'tGetCorrectedIC
  | Don'tGetCorrectedIC -> Fstream.nil
  | GiveBadVector n -> Fstream.of_list [Don'tGetCorrectedIC; GetCorrectedIC]
                       @@ Fstream.map (fun n -> GiveBadVector n)
                           (Fstream.filter ((<>) (Carray.length model.vec))
                              (shrink_nat n))

let shrink_cmd model = function
  | SolveNormal t ->
    begin
      match model.next_query_time with
      | Some t -> Fstream.singleton ({ model with last_query_time = t;
                                                  next_query_time = None; },
                                     SolveNormal t)
      | None -> Fstream.map
                  (fun t -> ({ model with last_query_time = t; },
                             SolveNormal t))
                  (shrink_query_time model.last_query_time t)
    end
  | GetRootInfo -> Fstream.nil
  | GetNRoots -> Fstream.nil
  | SetAllRootDirections dir ->
    Fstream.map
      (fun dir -> (model, SetAllRootDirections dir))
      (shrink_root_direction dir)
  | CalcIC_Y (t, ic_buf) ->
    assert (not (model.solving));
    let model = { model with solving = true } in
    Fstream.map (fun ic_buf -> (model, CalcIC_Y (t, ic_buf)))
      (shrink_ic_buf model ic_buf)
    @@ Fstream.map (fun t -> ({ model with next_query_time = Some t },
                              CalcIC_Y (t, ic_buf)))
        (shrink_query_time (model.last_query_time +. discrete_unit) t)
  | SetRootDirection dirs ->
    Fstream.map
      (fun dirs -> (model, SetRootDirection dirs))
      (shrink_array shrink_root_direction ~shrink_size:false dirs)

(* Fix up a command to be executable in a given state (= the model parameter).
   The incoming model may have fewer equations than when the command was
   generated, or have a different starting time.  *)
let fixup_cmd model = function
  | SolveNormal t ->
    begin
      match model.next_query_time with
      | Some t -> ({ model with last_query_time = t; next_query_time = None; },
                   SolveNormal t)
      | None -> let t = max t model.last_query_time in
                ({ model with last_query_time = t }, SolveNormal t)
    end
  | CalcIC_Y (t, GiveBadVector n) when n = Carray.length model.vec ->
    (* FIXME: is it OK to perform CalcIC_Y multiple times?  What about with an
       intervening solve?  Without an intervening solve?  *)
    ({ model with next_query_time = Some t },
     CalcIC_Y (t, GiveBadVector (if n = 0 then 1 else (n-1))))
  | GetRootInfo | GetNRoots | CalcIC_Y _ | SetRootDirection _
  | SetAllRootDirections _ as cmd -> (model, cmd)

let gen_cmds model =
  gen_1pass_list gen_cmd model

let shrink_cmds model =
  shrink_1pass_list shrink_cmd fixup_cmd model

(* Note that although the interpreter mutates the model record, the entry point
   of the interpreter creates a deep copy each time, so it's OK for the shrunk
   model to share structures with each other and with the original model.  *)
let shrink_model model cmds =
  let fixup_cmds new_model cmds =
    (new_model, fixup_list fixup_cmd new_model cmds)
  in
  let update_roots model cmds (i, roots) =
    let n = Array.length roots in
    let model = { model with
                  roots = roots;
                  root_info = if i < 0 then model.root_info
                              else Roots.create n;
                  root_dirs = if i < 0 then model.root_dirs
                              else Array.make n RootDirs.IncreasingOrDecreasing
                }
    in
    (* Fix up commands to cope with a shrunk root function set.  The set of
       equations is not shrunk.  *)
    let fixup_cmd_with_drop model cmd =
      match cmd with
      | SolveNormal _ | GetRootInfo | GetNRoots | CalcIC_Y _
      | SetAllRootDirections _ -> fixup_cmd model cmd
      | SetRootDirection root_dirs ->
        if Array.length root_dirs = 0
        then fixup_cmd model (SetRootDirection root_dirs)
        else
          let i = if i < Array.length root_dirs then i else 0 in
          fixup_cmd model (SetRootDirection (array_drop_elem root_dirs i))
    in
    (model,
     fixup_list (if i < 0 then fixup_cmd else fixup_cmd_with_drop) model cmds)
  and update_time model cmds t =
    let model = { model with last_query_time = t; last_tret = t; t0 = t; } in
    fixup_cmds model cmds
  in
  Fstream.map (update_time model cmds) (shrink_t0 model.t0)
  @@ Fstream.map
      (update_roots model cmds)
      (shorten_shrink_array
         (shrink_pair
            (shrink_root_time model.t0)
            (shrink_choice [| Roots.Rising; Roots.Falling |]))
         model.roots)

(* Scripts (model + command list) *)

let gen_script () =
  let model = gen_model () in
  let cmds  = gen_cmds model ()
  in (model, cmds)

(* As noted above, it's OK to not copy the model structure wholesale, but to
   copy only the parts that are shrunk or need to be fixed up.  *)
let shrink_neqs model cmds =
  (* Reduce commands that are dependent on number of equations.  *)
  let drop_from_cmd i = function
    | CalcIC_Y (t, GiveBadVector n) when n > 0 ->
      CalcIC_Y (t, GiveBadVector (n-1))
    | cmd -> cmd
  in
  let copy_resfn_drop i = function
    | ResFnLinear slopes -> ResFnLinear (carray_drop_elem slopes i)
    | ResFnExpDecay coefs -> ResFnExpDecay (carray_drop_elem coefs i)
  in
  let drop_eq i =
    let model =
      { model with
        resfn = copy_resfn_drop i model.resfn;
        vec = carray_drop_elem model.vec i;
        vec' = carray_drop_elem model.vec' i;
        t0 = model.t0;
        vec0 = carray_drop_elem model.vec0 i;
        vec'0 = carray_drop_elem model.vec'0 i;
      }
    in (model, List.map (drop_from_cmd i) cmds)
  in
  let n = Carray.length model.vec in
  Fstream.guard (n > 1) (Fstream.map drop_eq (Fstream.enum_then (n-1) (n-2) 0))

let shrink_script (model, cmds) =
  shrink_neqs model cmds
  @@ shrink_model model cmds
  @@ Fstream.map (fun cmds -> (model, cmds)) (shrink_cmds model cmds)

(* Pretty-printing and result comparison.  *)

let cmp_eps = ref 1e-5

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

let pp_results, dump_results, show_results, display_results,
  print_results, prerr_results =
  printers_of_pp (fun ?(prec=0) fmt rs ->
    if !read_write_invariance then pp_vlist pp_result fmt rs
    else
      (* List one result per line, with step numbers starting from 1.  However,
         the first result is the result of init and this will be tagged as
         such.  So the steps are numbered 1 to (length rs - 1). *)
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
          Format.fprintf fmt "Step %s: " (pad_show (i+1));
          pp_cmd fmt cmd)
           (Fstream.of_list cmds)))

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

let is_exn = function
  | Exn _ -> true
  | _ -> false
let not_exn x = not (is_exn x)

(* Compilation and execution *)

let copy_file ~from_file ~to_file () =
  let infile  = open_in from_file
  and outfile = open_out to_file in
  (try while true do output_char outfile (input_char infile) done
   with End_of_file -> ());
  close_in infile;
  close_out outfile

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
  | ResFnExpDecay coefs ->
    for i = 0 to Carray.length vec - 1 do
      vec.{i} <- vec0.{i} *. exp (-. coefs.{i} *. (t -. t0));
      vec'.{i} <- -. coefs.{i} *. vec.{i}
    done

(* Find the smallest root(s) in the range (t1..t2].  *)
let find_roots roots root_dirs t1 t2 =
  let n = Array.length roots in
  let pos i = fst roots.(i) in
  let valid i =
    match root_dirs.(i), snd roots.(i) with
    | RootDirs.IncreasingOrDecreasing, _ -> true
    | RootDirs.Increasing, Roots.Rising -> true
    | RootDirs.Increasing, Roots.Falling -> false
    | RootDirs.Decreasing, Roots.Rising -> false
    | RootDirs.Decreasing, Roots.Falling -> true
    | _, Roots.NoRoot -> assert false
  in
  if n = 0 then []
  else
    let record  = ref t2
    and holders = ref []
    in
    for i = 0 to n-1 do
      if t1 < pos i && pos i <= !record
      then
        begin
          if valid i then
            (if pos i < !record then
               (record := pos i;
                holders := [i])
             else
               holders := i::!holders)
        end
    done;
    !holders

(* Run a single command on the model.  *)
let model_cmd model = function
  | SolveNormal t ->
    (* NB: we don't model interpolation failures -- t will be monotonically
       increasing.  *)
    assert (model.last_query_time <= t);
    (* Undocumented behavior (sundials 2.5.0): solve_normal with t=t0 usually
       fails with "tout too close to t0 to start integration", but sometimes
       succeeds.  Whether it succeeds seems to be unpredictable.  *)
    if t = model.t0 then Any
    else
      let tret, flag =
        (* NB: roots that the solver is already sitting on are ignored, unless
           dt = 0.  *)
        let tstart =
          if t = model.last_tret
          then model.last_tret -. time_epsilon
          else model.last_tret
        in
        match find_roots model.roots model.root_dirs tstart t with
        | [] ->
          (* Undocumented behavior (sundials 2.5.0): a non-root return from
             solve_normal resets root info to undefined.  *)
          model.root_info_valid <- false;
          t, SolverResult Ida.Continue
        | (i::_) as is ->
          let tret = fst model.roots.(i) in
          model.root_info_valid <- true;
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
    if model.root_info_valid then RootInfo (Roots.copy model.root_info)
    else
      (* FIXME: this should be Exn Ida.IllInput or something like that.  *)
      Type (RootInfo (Roots.copy model.root_info))
  | GetNRoots ->
    Int (Array.length model.roots)
  | CalcIC_Y (_, GiveBadVector _) -> Exn (Invalid_argument "")
     (* Undocumented behavior (sundials 2.5.0): GetConsistentIC() and CalcIC()
        can only be called before IDASolve().  However, IDA only detects
        illicit calls of GetConsistentIC(), and for calls to CalcIC() after
        IDASolve(), it seems to initialize with garbage.

        The OCaml binding contains a flag to detect and prevent this case.

     *)
  | CalcIC_Y (_, _) when model.solving -> Exn Ida.IllInput
  | CalcIC_Y (_, ic_buf) ->
    (* The t carried in the command is just a hint which the exact solver
       doesn't need.  The time at which vec and vec' should be filled is
       model.last_tret.  *)
    begin
      exact_soln model.resfn model.t0        model.vec0 model.vec'0
                             model.last_tret model.vec  model.vec';
      (* Undocumented behavior (sundials 2.5.0): apparently it's OK to call
         CalcIC() on the same session without an intervening IDASolve() or
         IDAReInit().  Therefore, we don't set model.solving here.  *)
      match ic_buf with
      | Don'tGetCorrectedIC -> Unit (* The current time will be slightly
                                       ahead of t0.  By exactly how much
                                       is hard to say.  *)
      | GetCorrectedIC -> carray model.vec
      | GiveBadVector _ -> assert false
    end
  | SetAllRootDirections dir ->
    if Array.length model.roots = 0
    then Exn Ida.IllInput
    else (Array.fill model.root_dirs 0 (Array.length model.root_dirs) dir;
          Unit)
  | SetRootDirection dirs ->
    if Array.length model.roots = 0
    then Exn Ida.IllInput
    else
      (* The binding adjusts the root array length as needed.  *)
      let ndirs = Array.length dirs
      and nroots = Array.length model.roots in
      let dirs =
        if ndirs = nroots then dirs
        else
          let d = Array.make nroots RootDirs.IncreasingOrDecreasing in
          Array.blit dirs 0 d 0 (min nroots ndirs);
          d
      in
      (model.root_dirs <- dirs;
       Unit)

(* Report a summary of the initial model.  *)
let model_init model = Aggr [Float model.t0;
                             carray model.vec;
                             carray model.vec']

(* Run a list of commands on the model.  *)
let model_run (model, cmds) =
  let model = copy_model model in
  (* Evaluation order requires this let here.  *)
  let head = model_init model in
  Fstream.cons head (Fstream.map (model_cmd model) (Fstream.of_list cmds))


let verbose = ref false

let ida_test_case_driver model cmds =
  let just_cmp = ref false
  and marshal_results = ref false
  in
  Arg.parse
    [("--just-cmp", Arg.Set just_cmp,
      "don't print anything; just give the exit code");
     ("--marshal-results", Arg.Set marshal_results,
      "for internal use only");
     ("--read-write-invariance", Arg.Set read_write_invariance,
      "print data in a format that can be fed to ocaml toplevel")]
    (fun _ -> ()) "a test case generated by quickcheck";
  let step = ref 0 in
  let mismatches = ref [] in
  let expected_results = ref [] in
  let actual_results = ref [] in
  let do_cmd thunk =
    let expected = if !step = 0 then model_init model
                   else model_cmd model cmds.(!step - 1)
    and actual = try Lazy.force thunk with exn -> Exn exn in
    if !just_cmp then
      (if not (result_matches expected actual)
       then exit 1)
    else if !marshal_results
    then Marshal.to_channel stdout actual []
    else
      (actual_results := actual::!actual_results;
       expected_results := expected::!expected_results;
       if not (result_matches expected actual)
       then mismatches := (!step, expected, actual)::!mismatches);
    step := !step + 1
  in
  let finish () =
    if not !just_cmp && not !marshal_results then
      match !mismatches with
      | [] -> Printf.printf "OK, test successful.\n"; 0
      | ms ->
        let prerrf fmt = Format.fprintf Format.err_formatter fmt in
        let prerr_mismatch (i, exp, act) =
          prerrf "Result mismatch on step %d:@\n got@\n %s@\n \
                  but expected@\n %s@\n"
            i (show_result exp) (show_result act)
        in
        prerrf "Test failed.@\n[Reason]@\n";
        List.iter prerr_mismatch !mismatches;
        prerrf "@\n[Test Case]@\n";
        prerr_script (model, Array.to_list cmds);
        prerrf "@\n[Program Output]@\n";
        prerr_results (List.rev !actual_results);
        prerrf "@\n[Expected Output]@\n";
        prerr_results (List.rev !expected_results);
        1
    else 0
  in
  do_cmd, finish

let with_file_descr descr f =
  let ret = try f descr with exn -> Unix.close descr; raise exn in
  Unix.close descr;
  ret

(* Just checks if the generated code's self-test is successful.  *)
let prop_ida_ok ml_file_of_script script =
  compile (ml_file_of_script script);
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

let quickcheck_ida ml_file_of_script max_tests =
  let prop = prop_ida_ok ml_file_of_script in
  let err = Format.err_formatter in
  let fprintf = Format.fprintf err in
  let result =
    if !verbose then
      quickcheck gen_script shrink_script ~pp_input:pp_script prop max_tests
    else
      quickcheck gen_script shrink_script prop max_tests
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
    let exit_status = Unix.system !test_exec_file in
    assert (exit_status = Unix.WEXITED 1);
    Some script
  | Some (script, Falsified stat) ->
    (* Other exit statuses mean the test case crashed.  We'll run it again,
       this time letting it produce as much output as it can, and compare that
       to the expected results.  *)
    compile (ml_file_of_script script);
    copy_file ~from_file:(!test_exec_file ^ ".ml")
              ~to_file:!test_failed_file
              ();
    Format.fprintf err "Failed test code saved in %s.@\n@\n[Reason]@\nTest \
                        code crashed.@\n@\n[Test Case]@\n"
      !test_failed_file;
    pp_script err script;
    Format.fprintf err "@\n@\n[Program Output]@\n";
    Format.pp_print_flush err ();
    flush stderr;
    ignore (Unix.system
              (if !read_write_invariance
               then !test_exec_file ^ " --read-write-invariance"
               else !test_exec_file));
    flush stdout;
    flush stderr;
    Format.fprintf err "@\n[Expected Output]@\n";
    pp_results err (Fstream.to_list (model_run script));
    Format.pp_print_newline err ();
    Some script
