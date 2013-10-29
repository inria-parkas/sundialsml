(* Compile IDA test cases into OCaml that uses Ida_serial.  *)
module Ida = Ida_serial
module Carray = Ida_serial.Carray
module Roots = Ida.Roots
open Pprint
open Quickcheck
open Quickcheck_sundials
open Pprint_ida
open Quickcheck_ida
open Camlp4.PreCast
open Syntax
open Ast
open Expr_of
open Expr_of_sundials
open Expr_of_ida
open Expr_of_ida_model
module Camlp4aux = Camlp4aux.Make (Camlp4.PreCast.Syntax)
open Camlp4aux

let _loc = Loc.ghost

let semis when_empty ctor = function
  | [] -> when_empty
  | e::es -> ctor (List.fold_left (fun e1 e2 -> Ast.ExSem (_loc, e1, e2)) e es)

let expr_array es = semis <:expr<[||]>> (fun e -> Ast.ExArr (_loc, e)) es

let expr_list es =
  List.fold_right (fun e es -> <:expr<$e$::$es$>>) es <:expr<[]>>

let expr_seq es = semis <:expr<()>> (fun e -> Ast.ExSeq (_loc, e)) es

let rec expr_of_resfn_impl = function
  | ResFnLinear slopes ->
    (* forall i. vec'.{i} = slopes.{i} *)
    let neqs = Carray.length slopes in
    let go i = <:expr<res.{$`int:i$}
                         <- vec'.{$`int:i$} -. $`flo:slopes.{i}$>>
    in <:expr<fun t vec vec' res ->
               $expr_seq (List.map go (enum 0 (neqs-1)))$>>
  | ResFnExpDecay coefs ->
    (* forall i. vec'.{i} = - coefs.{i} * vec.{i} *)
    let neqs = Carray.length coefs in
    let go i = <:expr<res.{$`int:i$}
                         <- vec'.{$`int:i$}
                            +. $`flo:coefs.{i}$ *. vec.{$`int:i$}>>
    in <:expr<fun t vec vec' res ->
               $expr_seq (List.map go (enum 0 (neqs-1)))$>>
  | ResFnDie ->
    <:expr<fun _ _ _ _ ->
            failwith "exception raised on purpose from residual function"
            >>

let jac_expr_of_resfn get set neqs resfn =
  let rec go resfn =
    match resfn with
    | ResFnLinear slopes ->
      (* forall i. vec'.{i} - slopes.{i} = 0 *)
      fun i j -> if i = j then <:expr<$set$ jac ($`int:i$, $`int:j$) c>>
                 else <:expr<$set$ jac ($`int:i$, $`int:j$) 0.>>
    | ResFnExpDecay coefs ->
      (* forall i. vec'.{i} + coefs.{i} * vec.{i} = 0 *)
      fun i j -> if i = j then <:expr<$set$ jac ($`int:i$, $`int:j$)
                                                (c +. $`flo:coefs.{i}$)>>
                 else <:expr<$set$ jac ($`int:i$, $`int:j$) 0.>>
    | ResFnDie -> assert false
  and ixs = enum 0 (neqs - 1) in
  match resfn with
  | ResFnDie -> <:expr<fun _ -> failwith "exception raised on purpose from \
                                          dense jacobian function">>
  | ResFnLinear _ | ResFnExpDecay _ ->
     <:expr<fun jac_arg jac ->
            let c = jac_arg.Ida.jac_coef in
            $expr_seq
              (List.concat
                (List.map (fun i -> List.map (go resfn i) ixs) ixs))$>>

let set_jac model session =
  let neqs = Carray.length model.vec in
  match model.solver with
  | Ida.Dense -> let dense_get = <:expr<Ida.Densematrix.get>>
                 and dense_set = <:expr<Ida.Densematrix.set>>
                 in <:expr<Ida.Dls.set_dense_jac_fn $session$
                           $jac_expr_of_resfn dense_get dense_set neqs
                                              model.resfn$>>
  | Ida.Band _ | Ida.Sptfqmr _ | Ida.Spbcg _ | Ida.Spgmr _
  | Ida.LapackBand _ | Ida.LapackDense _ ->
    raise (Failure "linear solver not implemented")

let expr_of_roots model roots =
  let n = Array.length roots in
  let set i =
    match roots.(i) with
    | r, Roots.Rising -> <:expr<g.{$`int:i$} <- t -. $`flo:r$>>
    | r, Roots.Falling -> <:expr<g.{$`int:i$} <- $`flo:r$ -. t>>
    | _, Roots.NoRoot -> assert false
  in
  let f ss i = <:expr<$ss$; $set i$>> in
  if n = 0 then <:expr<Ida.no_roots>>
  else if model.root_fails then
    <:expr<($`int:n$,
            (fun _ _ _ _ ->
              failwith "exception raised on purpose from root function"))>>
  else
    <:expr<($`int:n$,
               (fun t vec vec' g ->
                  $Fstream.fold_left f (set 0)
                    (Fstream.enum 1 (n-1))$))>>

type primed = Primed | NotPrimed
let expr_of_calc_ic_impl fcn ic_bufs t =
  let go (badvecs, expr, ret) = function
    | Primed, GetCorrectedIC ->
      (badvecs, <:expr<$expr$ ~y':vec'>>, (<:expr<carray vec'>>::ret))
    | NotPrimed, GetCorrectedIC ->
      (badvecs, <:expr<$expr$ ~y:vec>>, (<:expr<carray vec>>::ret))
    | _, Don'tGetCorrectedIC -> (badvecs, expr, ret)
    | Primed, GiveBadVector n ->
      (and_binding badvecs <:binding<badvec' = Carray.create $`int:n$>>,
       <:expr<$expr$ ~y':badvec'>>,
       [])
    | NotPrimed, GiveBadVector n ->
      (and_binding badvecs <:binding<badvec = Carray.create $`int:n$>>,
       <:expr<$expr$ ~y:badvec>>,
       [])
  in
  match List.fold_left go (BiNil _loc, fcn, []) ic_bufs with
  | badvecs, fcn, [] -> smart_let _loc badvecs <:expr<$fcn$ $`flo:t$; Unit>>
  | badvecs, fcn, [r] -> smart_let _loc badvecs <:expr<$fcn$ $`flo:t$; $r$>>
  | badvecs, fcn, rs ->
    smart_let _loc badvecs <:expr<$fcn$ $`flo:t$;
                                  Aggr $expr_list (List.rev rs)$>>

(* Generate the test code that executes a given command.  *)
let expr_of_cmd_impl model = function
  | SolveNormalBadVector (t, n) ->
    <:expr<let tret, flag = Ida.solve_normal session $`flo:t$
                              (Carray.create $`int:n$)
                              (Carray.create $`int:n$) in
           Aggr [Float tret; SolverResult flag; carray vec; carray vec']>>
  | SolveNormal t ->
    <:expr<let tret, flag = Ida.solve_normal session $`flo:t$ vec vec' in
           Aggr [Float tret; SolverResult flag; carray vec; carray vec']>>
  | CalcIC_Y (t, ic_buf) ->
    expr_of_calc_ic_impl <:expr<Ida.calc_ic_y session>> [(NotPrimed, ic_buf)] t
  | CalcIC_YaYd' (t, ic_buf_y, ic_buf_y') ->
    smart_let _loc
      <:binding<vartypes = Ida.VarTypes.of_array
                             $expr_of_array expr_of_var_type
                                (vartypes_of_model model)$>>
      (expr_of_calc_ic_impl
         <:expr<Ida.calc_ic_ya_yd' session vartypes>>
         [(NotPrimed, ic_buf_y);
          (Primed, ic_buf_y')]
         t)
  | GetRootInfo ->
    <:expr<let roots = Ida.Roots.create (Ida.nroots session) in
           Ida.get_root_info session roots;
           RootInfo roots>>
  | GetNRoots ->
    <:expr<Int (Ida.nroots session)>>
  | SetAllRootDirections dir ->
    <:expr<Ida.set_all_root_directions session $expr_of_root_direction dir$;
           Unit>>
  | SetVarTypes ->
    <:expr<Ida.set_var_types session (Ida.VarTypes.of_array
           $expr_of_array expr_of_var_type (vartypes_of_model model)$);
           Unit>>
  | SetSuppressAlg b -> <:expr<Ida.set_suppress_alg session $`bool:b$; Unit>>
  | ReInit params ->
    let roots =
      match params.reinit_roots with
      | None -> Ast.ExNil _loc
      | Some r -> <:expr<~roots:$expr_of_roots model r$>>
    in
    let vec0 =
      match params.reinit_vec0_badlen with
      | None -> expr_of_carray params.reinit_vec0
      | Some n -> <:expr<Carray.create $`int:n$>>
    and vec'0 =
      match params.reinit_vec'0_badlen with
      | None -> expr_of_carray params.reinit_vec'0
      | Some n -> <:expr<Carray.create $`int:n$>>
    in
    <:expr<Ida.reinit session
           ~linsolv:$expr_of_linear_solver params.reinit_solver$
           $roots$
           $`flo:params.reinit_t0$
           $vec0$
           $vec'0$;
           Unit
           >>
  | SetRootDirection dirs ->
    <:expr<Ida.set_root_direction session
           $expr_array (List.map expr_of_root_direction (Array.to_list dirs))$;
           Unit>>
  | SetStopTime t ->
    <:expr<Ida.set_stop_time session $`flo:t$; Unit>>

let expr_of_cmds_impl model = function
  | [] -> <:expr<()>>
  | cmds ->
    let sandbox exp = <:expr<do_cmd (lazy $exp$)>> in
    expr_seq (List.map (fun cmd -> sandbox (expr_of_cmd_impl model cmd))
                cmds)

let _ =
  register_expr_of_exn (fun cont -> function
      | Invalid_argument msg -> <:expr<Invalid_argument $`str:msg$>>
      | Not_found -> <:expr<Not_found>>
      | Ida.IllInput -> <:expr<Ida.IllInput>>
      | Ida.TooMuchWork -> <:expr<Ida.TooMuchWork>>
      | Ida.TooClose -> <:expr<Ida.TooClose>>
      | Ida.TooMuchAccuracy -> <:expr<Ida.TooMuchAccuracy>>
      | Ida.ErrFailure -> <:expr<Ida.ErrFailure>>
      | Ida.ConvergenceFailure -> <:expr<Ida.ConvergenceFailure>>
      | Ida.LinearSetupFailure -> <:expr<Ida.LinearSetupFailure>>
      | Ida.LinearInitFailure -> <:expr<Ida.LinearInitFailure>>
      | Ida.LinearSolveFailure -> <:expr<Ida.LinearSolveFailure>>
      | Ida.FirstResFuncFailure -> <:expr<Ida.FirstResFuncFailure>>
      | Ida.RepeatedResFuncErr -> <:expr<Ida.RepeatedResFuncErr>>
      | Ida.UnrecoverableResFuncErr -> <:expr<Ida.UnrecoverableResFuncErr>>
      | Ida.RootFuncFailure -> <:expr<Ida.RootFuncFailure>>
      | Ida.BadK -> <:expr<Ida.BadK>>
      | Ida.BadT -> <:expr<Ida.BadT>>
      | Ida.BadDky -> <:expr<Ida.BadDky>>
      | exn -> cont exn
    )

let ml_of_script (model, cmds) =
  <:str_item<
    module Ida = Ida_serial
    module Carray = Ida.Carray
    open Quickcheck_sundials
    open Quickcheck_ida
    open Pprint
    let model = $expr_of_model model$
    let cmds = $expr_of_array expr_of_cmd (Array.of_list cmds)$
    let do_cmd, finish, err_handler = test_case_driver model cmds
    let _ =
      let vec  = $expr_of_carray model.vec0$
      and vec' = $expr_of_carray model.vec'0$ in
      let session = Ida.init
                  $expr_of_linear_solver model.solver$
                  $expr_of_resfn_impl model.resfn$
                  ~roots:$expr_of_roots model model.roots$
                  ~t0:$`flo:model.t0$
                  vec vec'
      in
      Ida.ss_tolerances session 1e-9 1e-9;
      $set_jac model <:expr<session>>$;
      Ida.set_err_handler_fn session err_handler;
      do_cmd (lazy (Aggr [Float (Ida.get_current_time session);
                          carray vec; carray vec']));
      $expr_of_cmds_impl model cmds$;
      exit (finish ())
   >>

let randseed =
  Random.self_init ();
  ref (Random.int ((1 lsl 30) - 1))

let ml_file_of_script script src_file =
  Camlp4.PreCast.Printers.OCaml.print_implem ~output_file:src_file
    (ml_of_script script);
  let chan = open_out_gen [Open_text; Open_append; Open_wronly] 0 src_file in
  Printf.fprintf chan "\n(* generated with random seed %d, test case %d *)\n"
    !randseed !test_case_number;
  close_out chan

;;
let _ =
  let max_tests = ref 50 in
  let options = [("--exec-file", Arg.Set_string test_exec_file,
                  "test executable name \
                   (must be absolute, prefixed with ./, or on path)");
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
    "randomly generate programs using IDA and check if they work as expected";

  Printf.printf "random generator seed value = %d\n" !randseed;
  flush stdout;
  Random.init !randseed;
  size := 1;
  quickcheck_script ml_file_of_script !max_tests

