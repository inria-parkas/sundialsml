module Carray = Sundials.Carray
open Pprint
open Pprint_ida
open Quickcheck
open Quickcheck_sundials
module Roots = Sundials.Roots
module RootDirs = Sundials.RootDirs
open Camlp4.PreCast.Syntax
open Expr_of

(* A state-machine model for IDA sessions.  *)
type model =
  {
    resfn : resfn_type;
    mutable solver : model_linear_solver;

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
    mutable root_fails : bool;          (* root function should fail *)
    vec : Carray.t;
    vec' : Carray.t;
    mutable t0 : float;
    vec0 : Carray.t;
    vec'0 : Carray.t;

    mutable stop_time : float;          (* infinity when unset *)

    (* Set during command generation to force the next query time to a
       particular value, if there is a query.  *)
    mutable next_query_time : float option;
  }
and model_linear_solver =
  (* Just like linear_solver, but avoids components that depend on the
     nvector type.  Optional callbacks are replaced by booleans to indicate
     whether they should be specified.  *)
  | MDense of bool
  | MLapackDense of bool
  | MBand of model_band_init
  | MLapackBand of model_band_init
  | MSpgmr of model_spils_init
  | MSpbcg of model_spils_init
  | MSptfqmr of model_spils_init
and model_band_init = { mmupper : int;
                        mmlower : int;
                        mband_jac : bool; }
and model_spils_init = { mmaxl : int;
                         mprec_setup_fn : bool;
                         mjac_times_vec_fn : bool; }
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
        | SetStopTime of float
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
    reinit_root_fails : bool;
    reinit_solver : model_linear_solver;
    reinit_vec0 : Carray.t;
    reinit_vec0_badlen : int option;
    reinit_vec'0 : Carray.t;
    reinit_vec'0_badlen : int option;
  }
and roots_spec = (float * Ida.Roots.root_event) array
and ic_buf = GetCorrectedIC
           | Don'tGetCorrectedIC
           | GiveBadVector of int (* Pass in a vector of incorrect size. *)
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
                         Ida.solver_result = solver_result
                        )
                 ~optional:(solving, consistent, next_query_time,
                            last_tret, root_info_valid)
         (* expr_of is derived in expr_of_ida_model.ml to avoid linking camlp4
            into test cases.  If you change the ~alias list above, you probably
            have to update the Makefile's expr_of_ida_model.ml target too.  *)
         )

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
    root_fails = m.root_fails;
    consistent = m.consistent;
    vec = Carray.of_carray m.vec;
    vec' = Carray.of_carray m.vec';
    t0 = m.t0;
    vec0 = Carray.of_carray m.vec0;
    vec'0 = Carray.of_carray m.vec'0;
    next_query_time = m.next_query_time;
    stop_time = m.stop_time;
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

(* Destructively update the model according to a command.  If the expected
   return value is Exn, this function will actually raise that exception.  The
   wrapper model_cmd converts it to Exn.  Keep in mind that for most exceptions
   the error message is not checked, apart from a few specific messages
   expected to be stable across versions, like "index out of bounds".  *)
let model_cmd_internal model = function
  | SolveNormalBadVector _ -> invalid_arg ""
  | SolveNormal query_time ->
    (* NB: we don't model interpolation failures -- t will be monotonically
       increasing.  *)
    assert (model.last_query_time <= query_time);
    model.last_query_time <- query_time;
    (* Behavior documented in source code (sundials 2.5.0): exception raised if
       stop time is before the first query time.  NB, this only checks the
       *first* query time.  *)
    if not model.solving && model.stop_time < model.last_tret
    then raise Ida.IllInput;
    let t =
      (* Account for stop time.  *)
      min query_time model.stop_time
    in
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
          t,
          SolverResult (if t = model.stop_time then Ida.StopTimeReached
                        else Ida.Continue),
          []
        | (i::_) as is ->
          let tret = fst model.roots.(i) in
          (tret,
           (* If the return time coincides with a root, then it's reasonable
              for the solver to prioritize either.  In real code, this is most
              likely determined by the state of the floating point error.  *)
           (if tret = t
            then Type (SolverResult Ida.Continue)
            else SolverResult Ida.RootsFound),
           is)
      in
      exact_soln model.resfn model.t0 model.vec0 model.vec'0
                                 tret model.vec  model.vec';
      (* The root and query time have to be updated here, after we've
         checked that the residual function doesn't throw an exception.  *)
      if Array.length model.roots > 0 && model.root_fails then
        failwith "exception raised on purpose from root function";
      (* Undocumented behavior (sundials 2.5.0): a non-root return from
         solve_normal sets root info to undefined.  *)
      model.root_info_valid <- (roots_to_update <> []);
      Roots.reset model.root_info;
      List.iter
        (fun i -> Roots.set model.root_info i (snd model.roots.(i)))
        roots_to_update;
      model.last_tret <- tret;
      model.solving <- true;
      (* Undocumented behavior (sundials 2.5.0): stop time is reset when it's
         hit.  *)
      if tret = model.stop_time then model.stop_time <- infinity;
      Aggr [Float tret; flag; carray model.vec; carray model.vec']
  | SetStopTime t ->
    (* Behavior documented in source code (sundials 2.5.0): can't set stop time
       in the past, if solving.  Otherwise the check is deferred until solution
       is polled.  *)
    if model.solving && t < model.last_tret
    then raise Ida.IllInput
    else (model.stop_time <- t; Unit)
  | GetRootInfo ->
    if model.root_info_valid then RootInfo (Roots.copy model.root_info)
    else
      (* FIXME: this should be Exn Ida.IllInput or something like that.  *)
      Type (RootInfo (Roots.copy model.root_info))
  | GetNRoots ->
    Int (Array.length model.roots)
  | CalcIC_Y (_, GiveBadVector _) -> invalid_arg ""
     (* Undocumented behavior (sundials 2.5.0): GetConsistentIC() and CalcIC()
        can only be called before IDASolve().  However, IDA only detects
        illicit calls of GetConsistentIC(), and for calls to CalcIC() after
        IDASolve(), it seems to initialize with garbage.

        The OCaml binding contains a flag to detect and prevent this case.

     *)
  | CalcIC_Y (_, _) when model.solving -> invalid_arg ""
  | CalcIC_Y (tout1, ic_buf) ->
    begin
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
    end
  | CalcIC_YaYd' (tout1, GiveBadVector _, _)
  | CalcIC_YaYd' (tout1, _, GiveBadVector _) -> Exn (Invalid_argument "")
  | CalcIC_YaYd' _ when model.solving ->
    (* Right now, the binding sets var types before invoking IDACalcIC().  *)
    model.vartypes_set <- true;
    invalid_arg ""
  | CalcIC_YaYd' (tout1, ic_buf_y, ic_buf_y') ->
    begin
      (* var types are set before calling CalcIC *)
      model.vartypes_set <- true;
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
    end
  | SetVarTypes -> model.vartypes_set <- true; Unit
  | SetSuppressAlg _ -> if model.vartypes_set then Unit else invalid_arg ""
  | SetAllRootDirections dir ->
    if Array.length model.roots = 0 then raise Ida.IllInput;
    Array.fill model.root_dirs 0 (Array.length model.root_dirs) dir;
    Unit
  | SetRootDirection dirs ->
    if Array.length model.roots = 0 then raise Ida.IllInput;
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
    model.root_dirs <- dirs;
    Unit
  | ReInit params ->
    if params.reinit_vec0_badlen <> None || params.reinit_vec'0_badlen <> None
    then invalid_arg "reinit: incorrect length"
    ;
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
       (* Undocumented behavior (sundials 2.5.0): IDARootInit() with the same
          number of roots as before does not reset the root directions.  *)
       if n <> Array.length model.root_dirs then
         (model.root_dirs <- Array.make n RootDirs.IncreasingOrDecreasing;
          model.root_info <- Roots.create n))
    ;
    model.root_info_valid <- false;
    Carray.blit params.reinit_vec0 model.vec;
    Carray.blit params.reinit_vec0 model.vec0;
    Carray.blit params.reinit_vec'0 model.vec';
    Carray.blit params.reinit_vec'0 model.vec'0;
    model.t0 <- params.reinit_t0;
    model.next_query_time <- None;
    Unit

let model_cmd model cmd =
  try model_cmd_internal model cmd
  with exn -> Exn exn

(* Non-destructively update a model by a command, and return the new model.  *)
let update_model model cmd =
  let model' = copy_model model in
  let _ = model_cmd model' cmd in
  model'

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

let gen_solver neqs =
  match Random.int 1 with
  | 0 -> if Sundials.blas_lapack_supported && Random.bool ()
         then MLapackDense (Random.bool ())
         else MDense (Random.bool ())
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
    solver = gen_solver neqs;
    solving = false;
    last_query_time = t0;
    last_tret = t0;
    consistent = true;
    vartypes_set = false;
    roots = roots;
    root_info = Roots.create num_roots;
    root_info_valid = false;
    root_dirs = Array.make num_roots RootDirs.IncreasingOrDecreasing;
    root_fails = Random.int 100 < 80;
    vec   = Carray.of_carray vec0;
    vec'  = Carray.of_carray vec'0;
    t0 = t0;
    vec0  = vec0;
    vec'0 = vec'0;
    next_query_time = None;
    stop_time = infinity;
  }

let shrink_ic_buf model = function
  | GetCorrectedIC -> Fstream.singleton Don'tGetCorrectedIC
  | Don'tGetCorrectedIC -> Fstream.nil
  | GiveBadVector n -> Fstream.of_list [Don'tGetCorrectedIC; GetCorrectedIC]
                       @+ Fstream.map (fun n -> GiveBadVector n)
                           (Fstream.filter ((<>) (Carray.length model.vec))
                              (shrink_nat n))

let shrink_solve_time model t =
  match model.next_query_time with
  | Some t' -> if t' < t
               then Fstream.singleton t'
               else Fstream.nil
  | None -> shrink_query_time (model.last_query_time +. discrete_unit) t

let fixup_ic_buf new_model_vec_len ic_buf =
  match ic_buf with
  | GiveBadVector n ->
    (* Ensure n doesn't match length of model's vector.  *)
    if n = new_model_vec_len
    then GiveBadVector (if n = 0 then 1 else (n-1))
    else ic_buf
  | Don'tGetCorrectedIC | GetCorrectedIC -> ic_buf

(* Fix up a command to be executable in a given state (= the model parameter).
   This function must be able to handle all fallouts from perturbations done to
   the model in shrink_model below.  This function returns just an updated hint
   along with the new command; the model is updated in the wrapper fixup_cmd
   below.
 *)
let fixup_just_cmd hint model cmd =
  match cmd with
  | ReInit params ->
    let hint' =
      match params.reinit_roots with
      | None -> hint
      | Some _ -> { hint_root_drop = None }
    in
    (hint', cmd)
  | SolveNormal t ->
    begin
      match model.next_query_time with
      | Some t -> (hint, SolveNormal t)
      | None -> let t = max t model.last_query_time in
                (hint, SolveNormal t)
    end
  | SolveNormalBadVector (t, n) when n = Carray.length model.vec ->
    (hint, SolveNormalBadVector (t, if n = 0 then 1 else (n-1)))
  | CalcIC_Y (tout1, ic_buf) ->
    (* tout1 is the time of the next query which may or may not exist.  If
       last_query_time has been changed, we have to bump up this future query
       time.  last_query_time however does not change after CalcIC_Y.  *)
    (* Is it OK to perform CalcIC_Y multiple times?  What about with an
       intervening solve?  Without an intervening solve?  *)
    let tout1' = max tout1 (model.last_query_time +. discrete_unit) in
    (hint, CalcIC_Y (tout1', fixup_ic_buf (Carray.length model.vec) ic_buf))
  | CalcIC_YaYd' (tout1, ic_buf_y, ic_buf_y') ->
    let tout1' = max tout1 (model.last_query_time +. discrete_unit) in
    let fixup_ic_buf = fixup_ic_buf (Carray.length model.vec) in
    (hint,
     CalcIC_YaYd' (tout1', fixup_ic_buf ic_buf_y, fixup_ic_buf ic_buf_y'))
  | SetRootDirection root_dirs as cmd ->
    let model_roots = Array.length model.roots
    and cmd_roots = Array.length root_dirs in
    if model_roots = cmd_roots then (hint, cmd)
    else
      (match hint.hint_root_drop with
       | None when cmd_roots < model_roots ->
         (* pad root_dirs *)
         (hint,
          SetRootDirection
            (Array.init model_roots
               (fun i -> if i < cmd_roots then root_dirs.(i)
                         else RootDirs.IncreasingOrDecreasing)))
       | None ->
         (* curtail root_dirs *)
         (hint, SetRootDirection (Array.sub root_dirs 0 model_roots))
       | Some _ when Array.length root_dirs = 0 -> (hint, cmd)
       | Some i ->
         let idrop = if i < Array.length root_dirs then i else 0 in
         (hint, SetRootDirection (array_drop_elem root_dirs idrop)))
  | SetVarTypes -> (hint, cmd)
  (* These commands need no fixing up, AND doesn't update the model in any
     way that is relevant to shrinking.  *)
  | SetSuppressAlg _ | GetRootInfo | GetNRoots | SetStopTime _
  | SolveNormalBadVector _ | SetAllRootDirections _ -> (hint, cmd)

let fixup_cmd (hint, model) cmd =
  let hint, cmd = fixup_just_cmd hint model cmd in
  let model = update_model model cmd in
  (hint, model), cmd

(* Commands: note that which commands can be tested without trouble depends on
   the state of the model.  *)


let gen_cmd =
  let _ () =
    (* This is a no-op that causes the compiler to direct you here whenever you
       add a new command.  Once you finish adding a generator for that command,
       add that command to this match.  *)
    match failwith "oops" with
    | SolveNormal _ | ReInit _ | SolveNormalBadVector _
    | CalcIC_Y _ | CalcIC_YaYd' _ | GetRootInfo | GetNRoots | SetVarTypes
    | SetSuppressAlg _ | SetAllRootDirections _ | SetRootDirection _
    | SetStopTime _
    -> ()
  in
  let always_predictable _ = true
  and gen_solve_normal model =
    let t =
      match model.next_query_time with
      | None -> gen_query_time model.last_query_time ()
      | Some t -> t
    in SolveNormal t
  and gen_reinit model =
    let neqs = Carray.length model.vec in
    let t0 = gen_t0 () in
    let vec0, vec'0  = init_vec_for neqs t0 model.resfn in
    let params =
      {
        reinit_t0 = gen_t0 ();
        reinit_roots = if Random.int 100 < 30 then None
                       else Some (gen_roots t0 ());
        reinit_root_fails = Random.int 100 < 30;
        reinit_solver = gen_solver neqs;
        reinit_vec0 = Carray.of_carray vec0;
        reinit_vec0_badlen = if Random.int 100 < 95 then None
                             else Some (gen_nat_avoiding neqs);
        reinit_vec'0 = Carray.of_carray vec'0;
        reinit_vec'0_badlen = if Random.int 100 < 95 then None
                              else Some (gen_nat_avoiding neqs);
      }
    in
    ReInit params
  in
  let calc_ic_ya_yd'_predictable model =
    match model.resfn with
    | ResFnLinear _ | ResFnExpDecay _ | ResFnDie -> true
  and gen_calc_ic_ya_yd' model =
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
    CalcIC_YaYd' (t, ic_buf_y, ic_buf_y')
  and calc_ic_y_predictable model =
    match model.resfn with
    | ResFnLinear _ -> false (* ic_calc_y doesn't work on ResFnLinear. *)
    | ResFnDie | ResFnExpDecay _ -> true
  and gen_calc_ic_y model =
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
    CalcIC_Y (t, ic_buf)
  and set_stop_time_predictable model =
    (* Undocumented behavior (sundials 2.5.0): set_stop_time fails if the stop
       time is already passed; however, the time is compared against IDA's
       internal time, not the last time that the user specified.  The internal
       time of IDA can't be guessed a priori, so if model.solving is true, we
       can't predict the outcome of set_stop_time.  *)
    not model.solving
  and gen_set_stop_time model = SetStopTime (gen_stop_time model.t0 ())
  in
  let choose =
    gen_cond_choice
      [| always_predictable, gen_solve_normal;
         always_predictable, gen_reinit;
         always_predictable, (fun model -> GetRootInfo);
         always_predictable, (fun model -> GetNRoots);
         always_predictable, (fun model -> SetVarTypes);
         set_stop_time_predictable, gen_set_stop_time;
         calc_ic_ya_yd'_predictable, gen_calc_ic_ya_yd';
         calc_ic_y_predictable, gen_calc_ic_y;
         always_predictable, (fun model -> SetSuppressAlg (gen_bool ()));

         always_predictable,
         (fun model ->
            match gen_solve_normal model with
            | SolveNormal t ->
              SolveNormalBadVector
                (t, gen_nat_avoiding (Carray.length model.vec))
            | _ -> assert false);

         always_predictable,
         (fun model ->
            (* 20% of the time, we (may) choose an incorrect size.  *)
            let size = if Random.int 100 < 20 then gen_nat ()
              else Array.length model.roots
            in
            let dirs = gen_array ~size:size gen_root_direction ()
            in SetRootDirection dirs);

         always_predictable,
         (fun model -> SetAllRootDirections (gen_root_direction ()));
      |]
  in
  fun model ->
    (* The application to the first model executes the gen_cond_choice, the
       application to the second model runs the selected generator.  *)
    let cmd = choose model model in
    (update_model model cmd, cmd)

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
  @+ Fstream.map (update_roots model cmds) (shrink_roots model.t0 model.roots)
  @+ Fstream.map (fun rf -> { model with root_fails = rf }, cmds)
      (Fstream.guard model.root_fails (Fstream.singleton false))

let gen_cmds model =
  gen_1pass_list gen_cmd model

(* Shrink a command so that its outcome is predictable when executed under the
   given model state.  Return a stream of pairs (hint, cmd') where cmd' is the
   shrunk command and hint is some hint for fixup_cmd about what has changed,
   if applicable.  The model sould NOT be updated; the wrapper shrink_cmd will
   do that.  *)
let shrink_just_cmd model = function
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
    (* TODO: shrink the solver.  *)
    (match params.reinit_roots with
     | None -> Fstream.nil
     | Some roots ->
       Fstream.cons ({ hint_root_drop = None },
                     ReInit { params with reinit_roots = None })
         (Fstream.map
            (fun (i, r) ->
               ({ hint_root_drop = if i < 0 then None else Some i },
                ReInit { params with reinit_roots = Some r }))
            (shrink_roots params.reinit_t0 roots)))
    @+
    Fstream.guard params.reinit_root_fails
      (Fstream.singleton
         (model_nohint, ReInit { params with reinit_root_fails = false }))
    @+
    Fstream.map
      (fun t0 -> (model_nohint, ReInit { params with reinit_t0 = t0 }))
      (shrink_t0 params.reinit_t0)
    @+
    (match params.reinit_vec0_badlen with
     | None -> Fstream.nil
     | Some l ->
       Fstream.map
         (fun l' -> (model_nohint,
                    ReInit { params with reinit_vec0_badlen = Some l' }))
         (shrink_nat_avoiding (Carray.length model.vec) l))
    @+
    (match params.reinit_vec'0_badlen with
     | None -> Fstream.nil
     | Some l ->
       Fstream.map
         (fun l' -> (model_nohint,
                    ReInit { params with reinit_vec'0_badlen = Some l' }))
         (shrink_nat_avoiding (Carray.length model.vec) l))
  | SolveNormalBadVector (t, n) ->
    Fstream.map (fun t -> (model_nohint, SolveNormalBadVector (t, n)))
      (shrink_solve_time model t)
    @+
    Fstream.map (fun n -> (model_nohint, SolveNormalBadVector (t, n)))
      (Fstream.filter ((<>) (Carray.length model.vec)) (shrink_nat n))
  | SolveNormal t ->
    Fstream.map (fun t -> (model_nohint, SolveNormal t))
      (shrink_solve_time model t)
  | SetAllRootDirections dir ->
    Fstream.map (fun dir -> (model_nohint, SetAllRootDirections dir))
      (shrink_root_direction dir)
  | CalcIC_Y (t, ic_buf) ->
    Fstream.map (fun ic_buf -> (model_nohint, CalcIC_Y (t, ic_buf)))
      (shrink_ic_buf model ic_buf)
    @+ Fstream.map (fun t -> (model_nohint, CalcIC_Y (t, ic_buf)))
        (shrink_query_time (model.last_query_time +. discrete_unit) t)
  | CalcIC_YaYd' (t, ic_buf_y, ic_buf_y') ->
    Fstream.map (fun ic_buf_y -> (model_nohint,
                                  CalcIC_YaYd' (t, ic_buf_y, ic_buf_y')))
      (shrink_ic_buf model ic_buf_y)
    @+ Fstream.map
        (fun ic_buf_y' -> (model_nohint,
                           CalcIC_YaYd' (t, ic_buf_y, ic_buf_y')))
        (shrink_ic_buf model ic_buf_y')
    @+ Fstream.map (fun t -> (model_nohint,
                              CalcIC_YaYd' (t, ic_buf_y, ic_buf_y')))
        (shrink_query_time (model.last_query_time +. discrete_unit) t)
  | SetRootDirection dirs ->
    Fstream.map
      (fun dirs -> (model_nohint, SetRootDirection dirs))
      (shrink_array shrink_root_direction ~shrink_size:false dirs)
  | SetStopTime t ->
    Fstream.map (fun t -> (model_nohint, SetStopTime t))
      (shrink_stop_time model.t0 t)
  | GetRootInfo | GetNRoots | SetVarTypes | SetSuppressAlg _ -> Fstream.nil

let shrink_cmd (hint, model) cmd =
  (* Shrink should be called on exactly one command, and all subsequent
     commands are handled by fixup_cmd.  *)
  assert (hint = model_nohint);
  Fstream.map (fun (hint, cmd) -> ((hint, update_model model cmd), cmd))
    (shrink_just_cmd model cmd)

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
  @+ shrink_model model cmds
  @+ Fstream.map (fun cmds -> (model, cmds)) (shrink_cmds model cmds)

let is_exn = function
  | Exn _ -> true
  | _ -> false
let not_exn x = not (is_exn x)

(* Report a summary of the initial model.  *)
let model_init model = Aggr [Float model.t0;
                             carray model.vec;
                             carray model.vec']

type ida_model = model
and ida_cmd = cmd

module IDALang =
  struct
    type model = ida_model
    and cmd = ida_cmd

    type script = model * cmd list

    let pp_model = pp_model
    let pp_cmd = pp_cmd

    let copy_model = copy_model

    let gen_script = gen_script
    let shrink_script = shrink_script

    let exec_cmd = model_cmd
    let model_init = model_init

    let result_matches = result_matches

    type ml_file_of_script = script -> string -> unit
  end

include TestImperativeLang (IDALang)

