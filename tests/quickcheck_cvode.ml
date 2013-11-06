module Carray = Sundials.Carray
open Pprint
open Pprint_cvode
open Quickcheck
open Quickcheck_sundials
module Roots = Sundials.Roots
module RootDirs = Sundials.RootDirs
open Camlp4.PreCast.Syntax
open Expr_of

(* A state-machine model for CVODE sessions.  *)
type model =
  {
    rhsfn : rhsfn_type;
    mutable lmm : Cvode.lmm;
    mutable iter : model_iter;

    (* Set when SolveNormal has been called on this model.  *)
    mutable solving : bool;

    mutable last_query_time : float;
    mutable last_tret : float;
    mutable roots : roots_spec;
    mutable root_dirs : Cvode.root_direction array;
    mutable root_info : Roots.t;
    mutable root_info_valid : bool;
    mutable root_fails : bool;          (* root function should fail *)
    vec : Carray.t;
    mutable t0 : float;
    vec0 : Carray.t;

    mutable stop_time : float;          (* infinity when unset *)
  }
and model_iter =
  | MFunctional
  | MNewton of model_linear_solver
and model_linear_solver =
  (* Just like linear_solver, but avoids components that depend on the nvector
     type.  Optional callbacks are replaced by booleans to indicate whether
     they should be specified.  *)
  | MDense of bool
  | MLapackDense of bool
  | MBand of model_bandrange * bool
  | MLapackBand of model_bandrange * bool
  | MDiag
and model_bandrange = { mmupper : int;
                        mmlower : int; }
and cmd = SolveNormal of float
        | SolveNormalBadVector of float * int
        | GetRootInfo
        | SetRootDirection of Cvode.root_direction array
        | SetAllRootDirections of Cvode.root_direction
        | GetNRoots
        | ReInit of reinit_params
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
    reinit_iter : model_iter;
    reinit_vec0 : Carray.t;
    reinit_vec0_badlen : int option;
  }
and roots_spec = (float * Cvode.Roots.root_event) array
and rhsfn_type =
  | RhsFnLinear of Carray.t (* The ODE is
                                 y'.{i} = coefs.{i}
                               with general solution
                                 y.{i} = coefs.{i} * t + C
                             *)
  | RhsFnExpDecay of Carray.t       (* The ODE is
                                         y'.{i} = coefs.{i} * y.{i}
                                       with general solution
                                         y.{i} = C exp (coefs.{i} * t) *)
  | RhsFnDie                 (* A rhsfn that always raises an exception.
                                For the purposes of finding the number of
                                variables, which ones are algebraic, and so on,
                                this behaves like the 1-variable system
                                  y' = 0.
                                Note a 0-variable system causes problems.  *)
deriving (pretty ~alias:(Carray.t = carray,
                         Cvode.Roots.root_event = root_event,
                         Cvode.root_direction = root_direction,
                         Roots.t = root_info,
                         Cvode.iter = iter,
                         Cvode.lmm = lmm,
                         Cvode.solver_result = solver_result)
                 ~optional:(last_tret, root_info_valid)
         (* expr_of is derived in expr_of_cvode_model.ml to avoid linking camlp4
            into test cases.  If you change the ~alias list above, you probably
            have to update the Makefile's expr_of_cvode_model.ml target too.  *)
         )

let carray x = Carray (Carray.of_carray x)

let model_nohint = { hint_root_drop = None }

let rec copy_rhsfn = function
  | RhsFnLinear slopes -> RhsFnLinear (Carray.of_carray slopes)
  | RhsFnExpDecay coefs -> RhsFnExpDecay (Carray.of_carray coefs)
  | RhsFnDie -> RhsFnDie
let copy_model m =
  {
    (* Everything is copied explicitly here instead of using { m with ... }
       because this way the compiler forces us to inspect every new field when
       it's added, so we won't forget to deep copy any of them.  *)
    rhsfn = copy_rhsfn m.rhsfn;
    lmm = m.lmm;
    iter = m.iter;
    solving = m.solving;
    last_query_time = m.last_query_time;
    last_tret = m.last_tret;
    roots = Array.copy m.roots;
    root_dirs = Array.copy m.root_dirs;
    root_info = Roots.copy m.root_info;
    root_info_valid = m.root_info_valid;
    root_fails = m.root_fails;
    vec = Carray.of_carray m.vec;
    t0 = m.t0;
    vec0 = Carray.of_carray m.vec0;
    stop_time = m.stop_time;
  }


(* Model interpretation *)

let rec exact_soln typ t0 vec0 t vec =
  match typ with
  | RhsFnLinear slopes ->
    for i = 0 to Carray.length vec - 1 do
      vec.{i} <- slopes.{i} *. (t -. t0);
    done
  | RhsFnExpDecay coefs ->
    for i = 0 to Carray.length vec - 1 do
      vec.{i} <- vec0.{i} *. exp (-. coefs.{i} *. (t -. t0));
    done
  | RhsFnDie ->
    raise (Failure "exception raised on purpose from rhs function")

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
    (* Bug (sundials 2.5.0): the comment in sundials source code says an
       exception is raised if stop time is before the first query time, but
       this check doesn't always succeed.  *)
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
          SolverResult (if t = model.stop_time then Cvode.StopTimeReached
                        else Cvode.Continue),
          []
        | (i::_) as is ->
          let tret = fst model.roots.(i) in
          (tret,
           (* If the return time coincides with a root, then it's reasonable
              for the solver to prioritize either.  In real code, this is most
              likely determined by the state of the floating point error.  *)
           (if tret = t
            then Type (SolverResult Cvode.Continue)
            else SolverResult Cvode.RootsFound),
           is)
      in
      exact_soln model.rhsfn model.t0 model.vec0
                                 tret model.vec;
      (* The root and query time have to be updated here, after we've
         checked that the rhs function doesn't throw an exception.  *)
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
      Aggr [Float tret; flag; Type (Carray (Carray.create 0))]
  | SetStopTime t ->
    (* Behavior documented in source code (sundials 2.5.0): can't set stop time
       in the past, if solving.  Otherwise the check is deferred until solution
       is polled.  *)
    if model.solving && t < model.last_tret
    then raise Cvode.IllInput
    else (model.stop_time <- t; Unit)
  | GetRootInfo ->
    if model.root_info_valid then RootInfo (Roots.copy model.root_info)
    else
      (* FIXME: this should be Exn Cvode.IllInput or something like that.  *)
      Type (RootInfo (Roots.copy model.root_info))
  | GetNRoots ->
    Int (Array.length model.roots)
  | SetAllRootDirections dir ->
    if Array.length model.roots = 0 then raise Cvode.IllInput;
    Array.fill model.root_dirs 0 (Array.length model.root_dirs) dir;
    Unit
  | SetRootDirection dirs ->
    if Array.length model.roots = 0 then raise Cvode.IllInput;
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
    if params.reinit_vec0_badlen <> None
    then invalid_arg "reinit: incorrect length"
    ;
    model.solving <- false;
    model.last_query_time <- params.reinit_t0;
    model.last_tret <- params.reinit_t0;
    (match params.reinit_roots with
     | None -> ()
     | Some r ->
       let n = Array.length r in
       model.roots <- Array.copy r;
       (* Undocumented behavior (sundials 2.5.0): CVODERootInit() with the same
          number of roots as before does not reset the root directions.  *)
       if n <> Array.length model.root_dirs then
         (model.root_dirs <- Array.make n RootDirs.IncreasingOrDecreasing;
          model.root_info <- Roots.create n))
    ;
    model.root_info_valid <- false;
    Carray.blit params.reinit_vec0 model.vec;
    Carray.blit params.reinit_vec0 model.vec0;
    model.t0 <- params.reinit_t0;
    Unit

let model_cmd model cmd =
  try model_cmd_internal model cmd
  with exn -> Exn exn

(* Non-destructively update a model by a command, and return the new model.  *)
let update_model model cmd =
  let model' = copy_model model in
  let _ = model_cmd model' cmd in
  model'

let gen_rhsfn t0 neqs () =
  let _ () =
    (* This is a no-op that causes the compiler to direct you here whenever you
       add a new kind of rhs function.  *)
    match RhsFnLinear (Carray.create 0) with
    | RhsFnLinear _ | RhsFnExpDecay _ | RhsFnDie -> ()
  in
  let safe_rhsfn =
    gen_choice
      [|
        (fun () ->
           let v = Carray.create neqs in
           for i = 0 to Carray.length v - 1 do
             v.{i} <- float_of_int (i+1)
           done;
           RhsFnLinear v);
        (fun () ->
           let v = Carray.create neqs in
           (* All coefficients must be non-negative.  Initially we had the i-th
              equation's coefficient be i+1, but this turned out to trigger
              Cvode.TooMuchWork, presumably because the higher-coefficient
              functions decay so fast that they underflow.  Setting them to all
              1's seems to let us avoid that issue.  *)
           for i = 0 to Carray.length v - 1 do
             v.{i} <- float_of_int 1
           done;
           RhsFnExpDecay v);
      |]
      ()
  in
  if Random.int 100 < 90
  then safe_rhsfn
  else RhsFnDie

let gen_iter neqs =
  match Random.int 3 with
  | 0 -> MFunctional
  | 1 -> if Sundials.blas_lapack_supported && Random.bool ()
         then MNewton (MLapackDense (Random.bool ()))
         else MNewton (MDense (Random.bool ()))
  | 2 -> MNewton MDiag
  | _ -> assert false

let gen_lmm neqs = gen_choice [| Cvode.Adams; Cvode.BDF |]

let gen_roots t0 () =
  let gen_root () =
    (gen_root_time t0 (),
     gen_choice [| Roots.Rising; Roots.Falling |])
  in
  match Random.int 3 with
  | 0 -> [||]
  | _ -> uniq_array (gen_array gen_root ())

(* Returns some initial vec0.  *)
let init_vec_for neqs t0 = function
  | RhsFnLinear slopes -> Carray.init neqs 0.
  | RhsFnExpDecay coefs -> Carray.init neqs 1.
  | RhsFnDie -> Carray.of_array [|0.|]

let gen_model () =
  let neqs  = min 10 (gen_pos ()) in
  let t0 = gen_t0 () in
  let rhsfn = gen_rhsfn t0 neqs () in
  let vec0  = init_vec_for neqs t0 rhsfn in
  let neqs  = Carray.length vec0 in
  let roots = gen_roots t0 () in
  let num_roots = Array.length roots in
  {
    rhsfn = rhsfn;
    lmm = gen_lmm neqs;
    iter = gen_iter neqs;
    solving = false;
    last_query_time = t0;
    last_tret = t0;
    roots = roots;
    root_info = Roots.create num_roots;
    root_info_valid = false;
    root_dirs = Array.make num_roots RootDirs.IncreasingOrDecreasing;
    root_fails = Random.int 100 < 80;
    vec   = Carray.of_carray vec0;
    t0 = t0;
    vec0  = vec0;
    stop_time = infinity;
  }

let shrink_solve_time model t =
  shrink_query_time (model.last_query_time +. discrete_unit) t

(* Fix up a command to be executable in a given state (= the model parameter).
   This function must be able to handle all fallouts from perturbations done to
   the model in shrink_model below.  This function returns just an updated hint
   along with the new command; the model is updated in the wrapper fixup_cmd
   below.
 *)
let fixup_just_cmd hint model cmd =
  match cmd with
  | ReInit _ when model.stop_time <> infinity ->
    (* Undocumented behavior (sundials 2.5.0): doing a reinit when a stop time
       is set, but not reached, can lead to a session that produces garbage
       solutions when polled. *)
    (hint, GetNRoots)
  | ReInit params ->
    let hint' =
      match params.reinit_roots with
      | None -> hint
      | Some _ -> { hint_root_drop = None }
    in
    (hint', cmd)
  | SolveNormal t -> (hint, SolveNormal (max t model.last_query_time))
  | SolveNormalBadVector (t, n) when n = Carray.length model.vec ->
    (hint, SolveNormalBadVector (t, if n = 0 then 1 else (n-1)))
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
  | SetStopTime t ->
    (* Note t0 can grow as a result of shrinking if a ReInit command that
       resets t0 to a smaller value is removed.  *)
    if t <= model.t0 then (hint, SetStopTime (model.t0 +. stop_time_offs))
    else (hint, cmd)
  (* These commands need no fixing up, AND doesn't update the model in any
     way that is relevant to shrinking.  *)
  | GetRootInfo | GetNRoots
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
    | GetRootInfo | GetNRoots | SetAllRootDirections _ | SetRootDirection _
    | SetStopTime _
    -> ()
  in
  let always_predictable _ = true
  and gen_solve_normal model =
    let t = gen_query_time model.last_query_time () in
    SolveNormal t
  and reinit_predictable model =
    (* Undocumented behavior (sundials 2.5.0): doing a reinit when a stop time
       is set, but not reached, can lead to a session that produces garbage
       solutions when polled.  *)
    model.stop_time = infinity
  and gen_reinit model =
    let neqs = Carray.length model.vec in
    let t0 = gen_t0 () in
    let vec0  = init_vec_for neqs t0 model.rhsfn in
    let params =
      {
        reinit_t0 = gen_t0 ();
        reinit_roots = if Random.int 100 < 30 then None
                       else Some (gen_roots t0 ());
        reinit_root_fails = Random.int 100 < 30;
        reinit_iter = gen_iter neqs;
        reinit_vec0 = Carray.of_carray vec0;
        reinit_vec0_badlen = if Random.int 100 < 95 then None
                             else Some (gen_nat_avoiding neqs);
      }
    in
    ReInit params
  in
  let set_stop_time_predictable model =
    (* Undocumented behavior (sundials 2.5.0): set_stop_time fails if the stop
       time is already passed; however, the time is compared against CVODE's
       internal time, not the last time that the user specified.  The internal
       time of CVODE can't be guessed a priori, so if model.solving is true, we
       can't predict the outcome of set_stop_time.  *)
    not model.solving
  and gen_set_stop_time model =
    (* Since stop time inhibits reinit (see reinit_predictable), we must ensure
       the stop time expires promptly.  *)
    SetStopTime (min (model.t0 +. stop_time_offs +. 3. *. discrete_unit)
                     (gen_stop_time model.t0 ()))
  in
  let choose =
    gen_cond_choice
      [| always_predictable, gen_solve_normal;
         reinit_predictable, gen_reinit;
         always_predictable, (fun model -> GetRootInfo);
         always_predictable, (fun model -> GetNRoots);
         set_stop_time_predictable, gen_set_stop_time;

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
      (* If the stop time is unset, t0 can take on any value.  Otherwise,
         we need to shrink to some t0 > stop time.  *)
      (if model.stop_time = infinity then shrink_t0 params.reinit_t0
       else
         shrink_query_time
           (model.stop_time -. stop_time_offs +. discrete_unit)
           params.reinit_t0)
    @+
    (match params.reinit_vec0_badlen with
     | None -> Fstream.nil
     | Some l ->
       Fstream.map
         (fun l' -> (model_nohint,
                     ReInit { params with reinit_vec0_badlen = Some l' }))
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
  | SetRootDirection dirs ->
    Fstream.map
      (fun dirs -> (model_nohint, SetRootDirection dirs))
      (shrink_array shrink_root_direction ~shrink_size:false dirs)
  | SetStopTime t ->
    Fstream.map (fun t -> (model_nohint, SetStopTime t))
      (shrink_stop_time model.t0 t)
  | GetRootInfo | GetNRoots -> Fstream.nil

let shrink_cmd (hint, model) cmd =
  (* Shrink should be called on exactly one command, and all subsequent
     commands are handled by fixup_cmd.  *)
  assert (hint = model_nohint);
  Fstream.map (fun (hint, cmd) -> ((hint, update_model model cmd), cmd))
    (shrink_just_cmd model cmd)

let shrink_cmds model =
  shrink_1pass_list ~shrink_tails_first:true
    shrink_cmd fixup_cmd (model_nohint, model)

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
    | ReInit params ->
      ReInit { params with
               reinit_vec0  = carray_drop_elem params.reinit_vec0 i;
             }
    | cmd -> cmd
  in
  let rec copy_rhsfn_drop i = function
    | RhsFnLinear slopes -> RhsFnLinear (carray_drop_elem slopes i)
    | RhsFnExpDecay coefs -> RhsFnExpDecay (carray_drop_elem coefs i)
    | RhsFnDie -> assert false (* RhsFnDie is treated as 1 equation and
                                  can't be shrunk. *)
  in
  let drop_eq i =
    let model =
      { model with
        rhsfn = copy_rhsfn_drop i model.rhsfn;
        vec = carray_drop_elem model.vec i;
        t0 = model.t0;
        vec0 = carray_drop_elem model.vec0 i;
      }
    in (model, List.map (drop_from_cmd i) cmds)
  in
  let n = Carray.length model.vec in
  Fstream.guard (n > 1) (Fstream.map drop_eq (Fstream.enum_then (n-1) (n-2) 0))

let shrink_script (model, cmds) =
  shrink_neqs model cmds
  @+ Fstream.map (fun cmds -> (model, cmds)) (shrink_cmds model cmds)
  @+ shrink_model model cmds

let is_exn = function
  | Exn _ -> true
  | _ -> false
let not_exn x = not (is_exn x)

(* Report a summary of the initial model.  *)
let model_init model = Aggr [Float model.t0;
                             carray model.vec]

type cvode_model = model
and cvode_cmd = cmd

module CvodeLang =
  struct
    type model = cvode_model
    and cmd = cvode_cmd

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

include TestImperativeLang (CvodeLang)

