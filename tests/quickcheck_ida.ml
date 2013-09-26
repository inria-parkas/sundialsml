module Carray = Sundials.Carray
open Pprint
open Quickcheck
open Quickcheck_sundials
module Roots = Sundials.Roots
module RootDirs = Sundials.RootDirs
open Camlp4.PreCast.Syntax
open Expr_of

(* Derive pretty-printers for types imported from the Ida module.  It would be
   ideal if we could have the camlp4 extension extract type definitions from
   ../ida.ml and do everything automatically, but that's too much work.  For
   now, we'll make do by manually copying the definitions; the types aren't
   supposed to change often, and the compiler will detect bit rot in this
   case, so it should be OK.  *)
module PPIda = struct
  open Ida (* deriving needs access to type names like linear_solver *)
  open VarTypes (* ditto *)

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
  deriving external_types (pretty ~prefix:Ida, expr_of ~prefix:Ida)

  type var_type = Algebraic | Differential
  deriving external_types (expr_of ~prefix:Ida.VarTypes)
end
include PPIda                           (* Include the derived functions. *)

(* A state-machine model for IDA sessions.  *)
type model =
  {
    resfn : resfn_type;
    mutable solver : Ida.linear_solver;

    (* Set when SolveNormal or CalcIC_* has been called on this model.  *)
    mutable solving : bool;

    (* Set when the current state vector satisfies the DAE.  Not used right
       now.  *)
    mutable consistent : bool;
    mutable vartypes_set : bool;
    mutable last_query_time : float;
    mutable last_tret : float;
    mutable roots : roots_spec;
    mutable root_dirs : Ida.root_direction array;
    mutable root_info : Roots.t;
    mutable root_info_valid : bool;
    vec : Carray.t;
    vec' : Carray.t;
    mutable t0 : float;
    vec0 : Carray.t;
    vec'0 : Carray.t;

    (* Set during command generation to force the next query time to a
       particular value, if there is a query.  *)
    mutable next_query_time : float option;
  }
and cmd = SolveNormal of float
        | SolveNormalBadVector of float * int
        | GetRootInfo
        | SetRootDirection of Ida.root_direction array
        | SetAllRootDirections of Ida.root_direction
        | GetNRoots
        | CalcIC_Y of float * ic_buf
        | ReInit of reinit_params
        | SetVarTypes
        | CalcIC_YaYd' of float * ic_buf * ic_buf
        | SetSuppressAlg of bool
and script = model * cmd list
(* When the model is shrunk, the commands have to be fixed up so that we don't
   run commands whose outcomes are unpredictable.  In most cases the model
   gives enough context information to avoid pathologies, but in some cases it
   helps to have more hints; for instance, if the root set is shrunk, knowing
   which root was dropped helps to drop the corresponding entry from every
   root-related command.  This type contains those additional hints.  However,
   note that in some cases the hint cannot be supplied, so the fixup function
   must be prepared to do something sensible without any hints.  *)
and fixup_hint =
  { hint_root_drop : int option;
  }
and reinit_params =
  {
    reinit_t0 : float;
    reinit_roots : roots_spec option;
    reinit_solver : Ida.linear_solver;
    reinit_vec0 : Carray.t;
    reinit_vec'0 : Carray.t;
  }
and roots_spec = (float * Ida.Roots.root_event) array
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
  | ResFnDie                 (* A resfn that always raises an exception.
                                For the purposes of finding the number of
                                variables, which ones are algebraic, and so on,
                                this behaves like the 1-variable system
                                  y' = 0.
                                Note a 0-variable system causes problems.  *)
deriving (pretty ~alias:(Carray.t = carray,
                         Ida.Roots.root_event = root_event,
                         Ida.root_direction = root_direction,
                         Roots.t = root_info,
                         Ida.linear_solver = linear_solver,
                         Ida.solver_result = solver_result)
                 ~optional:(Int, Float, solving, consistent, next_query_time,
                            last_tret, root_info_valid)
         ,expr_of ~alias:(Carray.t = carray,
                          Ida.Roots.root_event = root_event,
                          Ida.root_direction = root_direction,
                          Roots.t = root_info,
                          Ida.linear_solver = linear_solver,
                          Ida.solver_result = solver_result))

let carray x = Carray (Carray.of_carray x)

let model_nohint = { hint_root_drop = None }

let rec copy_resfn = function
  | ResFnLinear slopes -> ResFnLinear (Carray.of_carray slopes)
  | ResFnExpDecay coefs -> ResFnExpDecay (Carray.of_carray coefs)
  | ResFnDie -> ResFnDie
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
    vartypes_set = m.vartypes_set;
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


(* Model interpretation *)

let rec exact_soln typ t0 vec0 vec'0 t vec vec' =
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
  | ResFnDie ->
    raise (Failure "exception raised on purpose from residual function")

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

(* Find the algebraic vs. differential classification for all variables as an
   ordinary array.  *)
let vartypes_of_model model =
  match model.resfn with
  | ResFnLinear _ | ResFnExpDecay _ ->
    Array.make (Carray.length model.vec) Ida.VarTypes.Differential
  | ResFnDie -> [|Ida.VarTypes.Differential|]

(* Destructively update the model according to a command.  *)
let model_cmd model = function
  | SolveNormalBadVector _ -> Exn Ida.IllInput
  | SolveNormal t ->
    (* NB: we don't model interpolation failures -- t will be monotonically
       increasing.  *)
    assert (model.last_query_time <= t);
    (* Undocumented behavior (sundials 2.5.0): solve_normal with t=t0 usually
       fails with "tout too close to t0 to start integration", but sometimes
       succeeds.  Whether it succeeds seems to be unpredictable.  *)
    if t = model.t0 then Any
    else
      let tret, flag, roots_to_update =
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
          t, SolverResult Ida.Continue, []
        | (i::_) as is ->
          let tret = fst model.roots.(i) in
          (tret,
           (* If the queried time coincides with a root, then it's reasonable
              for the solver to prioritize either.  In real code, this is most
              likely determined by the state of the floating point error.  *)
           (if tret = t then Type (SolverResult Ida.Continue)
            else SolverResult Ida.RootsFound),
           is)
      in
      (try
         exact_soln model.resfn model.t0 model.vec0 model.vec'0
                                    tret model.vec  model.vec';
         (* The root and query time have to be updated here, after we've
            checked that the residual function doesn't throw an exception.  *)
         model.root_info_valid <- (roots_to_update <> []);
         Roots.reset model.root_info;
         List.iter
           (fun i -> Roots.set model.root_info i (snd model.roots.(i)))
           roots_to_update;
         model.last_query_time <- t;
         model.last_tret <- tret;
         model.solving <- true;
         Aggr [Float tret; flag; carray model.vec; carray model.vec']
       with Failure "exception raised on purpose from residual function" as exn ->
         Exn exn)
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
  | CalcIC_Y (tout1, ic_buf) ->
    begin
      try
        (* tout1 is just a hint to IDA that exact_soln doesn't need, and the
           time at which vec and vec' should be filled is model.last_tret.
           However, we need to run the solver up to t because that might
           trigger an exception.  Whether the exception actually triggers is
           unpredictable, so if it triggers in the model, the expected output
           is Any.  *)
        exact_soln model.resfn model.t0  model.vec0 model.vec'0
                                  tout1  model.vec  model.vec';
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
      with
        Failure "exception raised on purpose from residual function" -> Any
    end
  | CalcIC_YaYd' (tout1, GiveBadVector _, _)
  | CalcIC_YaYd' (tout1, _, GiveBadVector _) -> Exn (Invalid_argument "")
  | CalcIC_YaYd' _ when model.solving ->
    (* Right now, the binding sets var types before invoking IDACalcIC().  *)
    model.vartypes_set <- true;
    Exn Ida.IllInput
  | CalcIC_YaYd' (tout1, ic_buf_y, ic_buf_y') ->
    begin
      (* var types are set before calling CalcIC *)
      model.vartypes_set <- true;
      try
        (* tout1 is just a hint to IDA that exact_soln doesn't need, and the
           time at which vec and vec' should be filled is model.last_tret.
           However, we need to run the solver up to t because that might
           trigger an exception.  Whether the exception actually triggers is
           unpredictable, so if it triggers in the model, the expected output
           is Any.  *)
        exact_soln model.resfn model.t0  model.vec0 model.vec'0
                               tout1     model.vec  model.vec';
        exact_soln model.resfn model.t0        model.vec0 model.vec'0
                               model.last_tret model.vec  model.vec';
        (* Undocumented behavior (sundials 2.5.0): apparently it's OK to call
           CalcIC() on the same session without an intervening IDASolve() or
           IDAReInit().  Therefore, we don't set model.solving here.  *)
        match ic_buf_y, ic_buf_y' with
        | GetCorrectedIC, Don'tGetCorrectedIC -> carray model.vec
        | Don'tGetCorrectedIC, GetCorrectedIC -> carray model.vec'
        | GetCorrectedIC, GetCorrectedIC -> Aggr [carray model.vec;
                                                  carray model.vec']
        | Don'tGetCorrectedIC, Don'tGetCorrectedIC -> Unit
        | _, _ -> assert false
      with
        Failure "exception raised on purpose from residual function" -> Any
    end
  | SetVarTypes -> model.vartypes_set <- true; Unit
  | SetSuppressAlg _ -> if model.vartypes_set then Unit else Exn Ida.IllInput
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
  | ReInit params ->
    model.solver <- params.reinit_solver;
    model.solving <- false;
    model.consistent <- true;
    model.last_query_time <- params.reinit_t0;
    model.last_tret <- params.reinit_t0;
    (match params.reinit_roots with
     | None -> ()
     | Some r ->
       let n = Array.length r in
       model.roots <- Array.copy r;
       model.root_dirs <- Array.make n RootDirs.IncreasingOrDecreasing;
       model.root_info <- Roots.create n);
    model.root_info_valid <- false;
    Carray.blit params.reinit_vec0 model.vec;
    Carray.blit params.reinit_vec0 model.vec0;
    Carray.blit params.reinit_vec'0 model.vec';
    Carray.blit params.reinit_vec'0 model.vec'0;
    model.t0 <- params.reinit_t0;
    model.next_query_time <- None;
    Unit

let gen_resfn t0 neqs () =
  let _ () =
    (* This is a no-op that causes the compiler to direct you here whenever you
       add a new kind of residual function.  *)
    match ResFnLinear (Carray.create 0) with
    | ResFnLinear _ | ResFnExpDecay _ | ResFnDie -> ()
  in
  let safe_resfn =
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
           (* All coefficients must be non-negative.  Initially we had the i-th
              equation's coefficient be i+1, but this turned out to trigger
              Ida.TooMuchWork, presumably because the higher-coefficient
              functions decay so fast that they underflow.  Setting them to all
              1's seems to let us avoid that issue.  *)
           for i = 0 to Carray.length v - 1 do
             v.{i} <- float_of_int 1
           done;
           ResFnExpDecay v);
      |]
      ()
  in
  if Random.int 100 < 90
  then safe_resfn
  else ResFnDie

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
let init_vec_for neqs t0 = function
  | ResFnLinear slopes -> Carray.init neqs 0., Carray.of_carray slopes
  | ResFnExpDecay coefs ->
    let vec0 = Carray.init neqs 1. in
    let vec'0 = Carray.of_carray coefs in
    Carray.mapi (fun i vec0_i -> -. coefs.{i} *. vec0.{i}) vec'0;
    vec0, vec'0
  | ResFnDie -> Carray.of_array [|0.|], Carray.of_array [|0.|]

let gen_model () =
  let neqs  = min 10 (gen_pos ()) in
  let t0 = gen_t0 () in
  let resfn = gen_resfn t0 neqs () in
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
    vartypes_set = false;
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

let shrink_ic_buf model = function
  | GetCorrectedIC -> Fstream.singleton Don'tGetCorrectedIC
  | Don'tGetCorrectedIC -> Fstream.nil
  | GiveBadVector n -> Fstream.of_list [Don'tGetCorrectedIC; GetCorrectedIC]
                       @@ Fstream.map (fun n -> GiveBadVector n)
                           (Fstream.filter ((<>) (Carray.length model.vec))
                              (shrink_nat n))

let shrink_solve_time model t =
  match model.next_query_time with
  | Some t' -> assert (t' > model.last_query_time);
               if t' < t
               then Fstream.singleton ({ model with last_query_time = t';
                                                    next_query_time = None; },
                                       t')
               else Fstream.nil
  | None -> Fstream.map
              (fun t -> ({ model with last_query_time = t; }, t))
              (shrink_query_time (model.last_query_time +. discrete_unit) t)

let fixup_ic_buf new_model_vec_len ic_buf =
  match ic_buf with
  | GiveBadVector n ->
    (* Ensure n doesn't match length of model's vector.  *)
    if n = new_model_vec_len
    then GiveBadVector (if n = 0 then 1 else (n-1))
    else ic_buf
  | Don'tGetCorrectedIC | GetCorrectedIC -> ic_buf

(* Fix up a command to be executable in a given state (= the model parameter).
   The incoming model will always have the same residual function as when the
   command was generated, but it may have a different starting time, or a
   different number of root functions.
   This function must be able to handle all fallouts from perturbations done to
   the model in shrink_model below.
   The model is updated only to the extent necessary to track which commands
   are testable.
 *)
let fixup_cmd ((diff, model) as ctx) cmd =
  match cmd with
  | ReInit params ->
    let model' = copy_model model
    and params' = { params with
                    reinit_solver = model.solver;
                  }
    and diff =
      match params.reinit_roots with
      | None -> diff
      | Some _ -> { hint_root_drop = None }
    in
    let cmd' = ReInit params' in
    ignore (model_cmd model' cmd');
    ((diff, model'), cmd')
  | SolveNormal t ->
    begin
      match model.next_query_time with
      | Some t -> ((diff, { model with last_query_time = t;
                                       next_query_time = None; }),
                   SolveNormal t)
      | None -> let t = max t model.last_query_time in
                ((diff, { model with last_query_time = t }), SolveNormal t)
    end
  | SolveNormalBadVector (t, n) when n = Carray.length model.vec ->
    ((diff, model), SolveNormalBadVector (t, if n = 0 then 1 else (n-1)))
  | SetSuppressAlg _ as cmd -> (ctx, cmd)
  | CalcIC_Y (t, ic_buf) ->
    (* Keep in mind the t carried in CalcIC_Y is the time of the next query
       which may or may not exist.  If last_query_time has been changed, we may
       have to bump up this future query time.  The last_query_time itself
       however does not change after CalcIC_Y.  *)
    (* Is it OK to perform CalcIC_Y multiple times?  What about with an
       intervening solve?  Without an intervening solve?  *)
    let t' = t (* max t (model.last_query_time +. discrete_unit) *) in
    ((diff, { model with next_query_time = Some t' }),
     CalcIC_Y (t', fixup_ic_buf (Carray.length model.vec) ic_buf))
  | CalcIC_YaYd' (t, ic_buf_y, ic_buf_y') ->
    let t' = max t (model.last_query_time +. discrete_unit) in
    let fixup_ic_buf = fixup_ic_buf (Carray.length model.vec) in
    ((diff, { model with next_query_time = Some t';
                         vartypes_set = true }),
     CalcIC_YaYd' (t', fixup_ic_buf ic_buf_y, fixup_ic_buf ic_buf_y'))
  | SetRootDirection root_dirs as cmd ->
    let model_roots = Array.length model.roots
    and cmd_roots = Array.length root_dirs in
    if model_roots = cmd_roots then (ctx, cmd)
    else
      (match diff.hint_root_drop with
       | None when cmd_roots < model_roots ->
         (* pad root_dirs *)
         (ctx,
          SetRootDirection
            (Array.init model_roots
               (fun i -> if i < cmd_roots then root_dirs.(i)
                         else RootDirs.IncreasingOrDecreasing)))
       | None ->
         (* curtail root_dirs *)
         (ctx,
          SetRootDirection
            (Array.sub root_dirs 0 model_roots))
       | Some _ when Array.length root_dirs = 0 -> (ctx, cmd)
       | Some i ->
         let idrop = if i < Array.length root_dirs then i else 0 in
         (ctx, SetRootDirection (array_drop_elem root_dirs idrop)))
  | GetRootInfo | GetNRoots | SetVarTypes
  | SolveNormalBadVector _ | SetAllRootDirections _ as cmd -> (ctx, cmd)

(* Commands: note that which commands can be tested without trouble depends on
   the state of the model.  *)

let gen_cmd =
  let _ () =
    (* This is a no-op that causes the compiler to direct you here whenever you
       add a new kind of command.  Once you finish adding a generator for a
       command, add that command's case to this match.  *)
    match failwith "oops" with
    | SolveNormal _ | ReInit _ | SolveNormalBadVector _
    | CalcIC_Y _ | CalcIC_YaYd' _ | GetRootInfo | GetNRoots | SetVarTypes
    | SetSuppressAlg _ | SetAllRootDirections _ | SetRootDirection _
    -> ()
  in
  let gen_solve_normal model =
    let t =
      match model.next_query_time with
      | None -> gen_query_time model.last_query_time ()
      | Some t -> t
    in ({ model with last_query_time = t }, SolveNormal t)
  and gen_reinit_params model =
    let neqs = Carray.length model.vec in
    let t0 = gen_t0 () in
    let vec0, vec'0  = init_vec_for neqs t0 model.resfn in
    let params =
      {
        reinit_t0 = gen_t0 ();
        reinit_roots = if Random.int 100 < 30 then None
                       else Some (gen_roots t0 ());
        reinit_solver = gen_solver false neqs;
        reinit_vec0 = Carray.of_carray vec0;
        reinit_vec'0 = Carray.of_carray vec'0;
      }
    in
    let model = copy_model model in
    ignore (model_cmd model (ReInit params));
    (model, ReInit params)
  and gen_nat_avoiding k =
    let n = gen_pos () in
    if n <= k then n-1
    else n
  in
  let gen_calc_ic_ya_yd' model =
    match model.resfn with
    | ResFnLinear _ | ResFnExpDecay _ | ResFnDie ->
      let t = gen_query_time (model.last_query_time +. discrete_unit) ()
      and gen_ic_buf_y, gen_ic_buf_y' =
        let b () = GiveBadVector
            (gen_nat_avoiding (Carray.length model.vec))
        and c () = GetCorrectedIC
        and d () = Don'tGetCorrectedIC in
        gen_weighted_choice
          [| (30, (c,c));  (10, (c,d));  (10, (c,b));
             (10, (d,c));  ( 5, (d,d));  (10, (d,b));
             (10, (b,c));  (10, (b,d));  ( 5, (b,b)); |]
          ()
      in
      let ic_buf_y  = gen_ic_buf_y ()
      and ic_buf_y' = gen_ic_buf_y' () in
      ({ model with next_query_time = Some t },
       CalcIC_YaYd' (t, ic_buf_y, ic_buf_y'))
  and gen_calc_ic_y model =
    match model.resfn with
    | ResFnLinear _ ->
      (* ic_calc_y doesn't work on ResFnLinear. *)
      (model, GetRootInfo)
    | ResFnDie | ResFnExpDecay _ ->
      let t = gen_query_time (model.last_query_time +. discrete_unit) ()
      and rand = Random.int 100 in
      let ic_buf =
        (* With 40% probability, don't receive corrected vector.
           With 10% probability, pass in bad vector.  *)
        if rand < 40 then Don'tGetCorrectedIC
        else if rand < 50
        then GiveBadVector (gen_nat_avoiding (Carray.length model.vec))
        else GetCorrectedIC
      in
      ({ model with next_query_time = Some t }, CalcIC_Y (t, ic_buf))
  in
  fun model ->
    gen_choice
      [| gen_solve_normal;
         gen_reinit_params;
         (fun model ->
            match gen_solve_normal model with
            | (model, SolveNormal t) ->
              (model, SolveNormalBadVector
                 (t, (gen_nat_avoiding (Carray.length model.vec))))
            | _ -> assert false);
         (fun model -> (model, GetRootInfo));
         (fun model -> (model, GetNRoots));
         (fun model -> (model, SetVarTypes));
         (fun model ->
            (* 20% of the time, we (may) choose an incorrect size.  *)
            let size = if Random.int 100 < 20 then gen_nat ()
              else Array.length model.roots
            in
            let dirs = gen_array ~size:size gen_root_direction ()
            in (model, SetRootDirection dirs));
         (fun model -> (model, SetAllRootDirections (gen_root_direction ())));
         gen_calc_ic_ya_yd';
         gen_calc_ic_y;
         (fun model -> model, SetSuppressAlg (gen_bool ()));
      |]
      model


let shrink_roots t0 =
  shorten_shrink_array
    (shrink_pair
       (shrink_root_time t0)
       (shrink_choice [| Roots.Rising; Roots.Falling |]))

(* Shrink a script by simplifying the model, but without removing equations.
   The set of changes attempted here must be kept in sync with the set of
   changes that fixup_cmd above can handle.

   Note that although the interpreter mutates the model record, the entry point
   of the interpreter creates a deep copy each time, so it's OK for the shrunk
   model to share structures with each other and with the original model.  *)
let shrink_model model cmds =
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
    let diff = if i < 0 then model_nohint else { hint_root_drop = Some i } in
    (model, fixup_list fixup_cmd (diff, model) cmds)
  and update_time model cmds t =
    let model = { model with last_query_time = t; last_tret = t; t0 = t; } in
    (model, fixup_list fixup_cmd (model_nohint, model) cmds)
  in
  Fstream.map (update_time model cmds) (shrink_t0 model.t0)
  @@ Fstream.map (update_roots model cmds) (shrink_roots model.t0 model.roots)

let gen_cmds model =
  gen_1pass_list gen_cmd model

let shrink_cmd ((diff, model) as ctx) cmd =
  (* Either shrink the model or shrink a command; shouldn't have to do both or
     shrink multiple commands at the same time.  *)
  assert (diff = model_nohint);
  match cmd with
  | ReInit params ->
    (* Let p = incoming params, p' = outgoing params,
           d = incoming diff,   d' = outgoing diff.
       Abbreviate hint_root_drop as root_drop, reinit_roots as roots, and Some
       x as x.  Then d'.root_drop is determined by the following table.  Keep
       in mind None can mean "no info" rather than "root set unchanged".

       p.roots      =    _    |       _        |     None    |   Some _
       p'.roots     = p.roots | drop i p.roots |     None    |   None
       -----------------------------------------------------------------
       d'.root_drop =   None  |     Some i     | d.root_drop |   None
    *)
    let ret diff' params' =
      let model' = copy_model model in
      let cmd = ReInit params' in
      ignore (model_cmd model' cmd);
      ((diff', model'), cmd)
    in
    (* TODO: implement solver shrinking.  *)
    (match params.reinit_roots with
     | None -> Fstream.nil
     | Some roots ->
       Fstream.cons (ret { hint_root_drop = None }
                         { params with reinit_roots = None })
         (Fstream.map
            (fun (i, r) ->
               ret { hint_root_drop = Some i }
                   { params with reinit_roots = Some r })
            (shrink_roots params.reinit_t0 roots)))
    @@
    Fstream.map (fun t0 -> ret diff { params with reinit_t0 = t0 })
      (shrink_t0 params.reinit_t0)
  | SolveNormalBadVector (t, n) ->
    Fstream.map (fun (m,t) -> ((diff, m), SolveNormalBadVector (t, n)))
      (shrink_solve_time model t)
    @@
    Fstream.map (fun n -> (ctx, SolveNormalBadVector (t, n)))
      (Fstream.filter ((<>) (Carray.length model.vec)) (shrink_nat n))
  | SolveNormal t ->
    Fstream.map (fun (m,t) -> ((diff, m), SolveNormal t))
      (shrink_solve_time model t)
  | SetAllRootDirections dir ->
    Fstream.map (fun dir -> (ctx, SetAllRootDirections dir))
      (shrink_root_direction dir)
  | CalcIC_Y (t, ic_buf) ->
    assert (not model.solving);
    let model' = { model with solving = true; next_query_time = Some t } in
    let ctx' = (diff, model') in
    Fstream.map (fun ic_buf -> (ctx', CalcIC_Y (t, ic_buf)))
      (shrink_ic_buf model ic_buf)
    @@ Fstream.map (fun t -> (ctx', CalcIC_Y (t, ic_buf)))
        (shrink_query_time (model.last_query_time +. discrete_unit) t)
  | CalcIC_YaYd' (t, ic_buf_y, ic_buf_y') ->
    assert (not model.solving);
    let model' = { model with solving = true; next_query_time = Some t } in
    let ctx' = (diff, model') in
    Fstream.map (fun ic_buf_y -> (ctx', CalcIC_YaYd' (t, ic_buf_y, ic_buf_y')))
      (shrink_ic_buf model ic_buf_y)
    @@ Fstream.map
        (fun ic_buf_y' -> (ctx', CalcIC_YaYd' (t, ic_buf_y, ic_buf_y')))
        (shrink_ic_buf model ic_buf_y')
    @@ Fstream.map (fun t -> (ctx', CalcIC_YaYd' (t, ic_buf_y, ic_buf_y')))
        (shrink_query_time (model.last_query_time +. discrete_unit) t)
  | SetRootDirection dirs ->
    Fstream.map
      (fun dirs -> (ctx, SetRootDirection dirs))
      (shrink_array shrink_root_direction ~shrink_size:false dirs)
  | GetRootInfo | GetNRoots | SetVarTypes | SetSuppressAlg _ -> Fstream.nil

let shrink_cmds model =
  shrink_1pass_list shrink_cmd fixup_cmd (model_nohint, model)

(* Scripts (model + command list) *)

let gen_script () =
  let model = gen_model () in
  let cmds  = gen_cmds model ()
  in (model, cmds)

(* Shrink a script by removing an equation from it.
   As noted above, it's OK to not copy the model structure wholesale, but to
   copy only the parts that are shrunk or need to be fixed up.  *)
let shrink_neqs model cmds =
  (* Reduce commands that are dependent on number of equations.  *)
  let drop_from_cmd i = function
    | CalcIC_Y (t, GiveBadVector n) when n > 0 ->
      CalcIC_Y (t, GiveBadVector (n-1))
    | ReInit params ->
      ReInit { params with
               reinit_vec0  = carray_drop_elem params.reinit_vec0 i;
               reinit_vec'0 = carray_drop_elem params.reinit_vec'0 i;
             }
    | cmd -> cmd
  in
  let rec copy_resfn_drop i = function
    | ResFnLinear slopes -> ResFnLinear (carray_drop_elem slopes i)
    | ResFnExpDecay coefs -> ResFnExpDecay (carray_drop_elem coefs i)
    | ResFnDie -> assert false (* ResFnDie is treated as 1 equation and
                                  can't be shrunk. *)
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
          prerrf "Result mismatch on step %d:@\ngot@\n  %s@\n\
                  but expected@\n  %s@\n"
            i (show_result act) (show_result exp)
        in
        prerrf "Test failed.@\n[Reason]@\n";
        List.iter prerr_mismatch !mismatches;
        prerrf "@\n[Test Case]@\n";
        prerr_script (model, Array.to_list cmds);
        prerrf "@\n@\n[Program Output]@\n";
        prerr_results (List.rev !actual_results);
        prerrf "@\n@\n[Expected Output]@\n";
        prerr_results (List.rev !expected_results);
        prerrf "@\n";
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
    Format.fprintf err "@\n[Reason]@\nTest \
                        code crashed.@\n@\n[Test Case]@\n";
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
    Format.fprintf err
      "@\n@\nFailed test code saved in %s, compile with:@\n%s@\n"
      !test_failed_file (!test_compiler ^ " " ^ !test_failed_file);
    Some script
