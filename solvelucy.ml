(*
 * Timothy Bourke (INRIA) & Marc Pouzet (ENS), August 2009
 *)

module Cvode = Cvode_serial
module Roots = Cvode.Roots
module Carray = Cvode.Carray

type lucyf =
   bool
  -> Roots.t
  -> Cvode.val_array
  -> Cvode.der_array
  -> Cvode.rootval_array
  -> bool

let lmm = ref Cvode.Adams
let iter = ref Cvode.Functional
let step = ref Cvode.normal

let show_root_names = ref false

let max_sim_time = ref None
let min_step_size = ref None
let max_step_size = ref (0.1)

let rel_tol = ref None
let abs_tol = ref None

let epsilons_to_add = ref 0

exception TooManyZeroCrossings

let zeroc_limit = ref 0
let zeroc_ignore_at_limit = ref false
let zeroc_deadzone = ref 0.0
let zeroc_time_tol = ref epsilon_float

let log = ref false
let log_zeroc = ref false

let set_zeroc_limit n = (zeroc_limit := n)
let ignore_at_zeroc_limit () = (zeroc_ignore_at_limit := true)
let set_zeroc_deadzone z = (zeroc_deadzone := z)

let enable_logging () = (log := true)
let enable_zeroc_logging () = (log_zeroc := true)

let printf = Printf.printf

let add_epsilons num_eps v =
  v +. (float(num_eps)
        *. (if v = 0.0 then min_float else abs_float v)
        *. epsilon_float
        *. 100.0)

let run allow_delta (lf : lucyf) advtime states roots =
  let n_roots = Array.length roots
  and n_cstates = Array.length states
  in

  let cstates    = Carray.create n_cstates
  and cder       = Carray.create n_cstates

  and roots_in   = Roots.create n_roots
  and no_roots_in = Roots.create n_roots

  and roots_out  = Carray.create n_roots
  and roots_out' = Carray.create n_roots

  and last_discrete_t = ref (-0.0)
  and discrete_tally = ref (0)

  and advtime =
    match advtime with
    | None -> (fun t -> t +. !max_step_size)
    | Some f -> f
  in

  let f t cs ds =
    ignore (lf false roots_in cs ds roots_out)
  and g t cs rs =
    ignore (lf false no_roots_in cs cder rs);
    Carray.clamp !zeroc_deadzone rs;

    if !log_zeroc then begin
      print_endline "(----";
      print_string " ZC:";
      Carray.print_with_time t cs;
      print_string " ZR:";
      Carray.print_with_time t rs;
      print_endline " ----)"
    end
  in

  let calculate_roots_out t = g t cstates roots_out in

  (* calculate ri by comparing ro (before) to ro (after). *)
  let calculate_roots_in ro ro' =
    let rin = Roots.set roots_in in
    for i = 0 to Carray.length ro - 1 do
      rin i (ro.{i} < 0.0 && ro'.{i} >= 0.0);
    done
  in

  let rec init () =
    Carray.fill cder 0.0;

    (* INIT CALL *)
    ignore (lf true roots_in cstates cder roots_out);
    if !epsilons_to_add <> 0 then
      Carray.map (add_epsilons !epsilons_to_add) cstates;

    let s = Cvode.init (!lmm) (!iter) f (n_roots, g) cstates in
    Cvode.set_all_root_directions s Cvode.Increasing;

    (match !max_sim_time with None -> () | Some t -> Cvode.set_stop_time s t);
    (match !min_step_size with None -> () | Some t -> Cvode.set_min_step s t);
    (match !rel_tol, !abs_tol with
     | Some rt, Some at -> Cvode.ss_tolerances s rt at | _ -> ());

    Roots.reset roots_in;
    if !log then begin
      print_string "H : time";
      Array.iter (printf "\t%s") states;
      print_newline ();

      print_string "--+\n";
      print_string "I : ";
      Carray.print_with_time 0.0 cstates
    end;
    continuous s (advtime 0.0)

  and continuous s t =
    (* CONTINUOUS CALL(S) *)
    (* INV: forall i. roots_in[i] = false *)
    let (t', result) = !step s t cstates
    in
      if !log then begin
        if result = Cvode.RootsFound
        then print_string "C': "
        else print_string "C : ";
        Carray.print_with_time t' cstates
      end;
      match result with
      | Cvode.RootsFound -> begin
            Cvode.get_root_info s roots_in;
            calculate_roots_out t';
            (* NB: we are forced to recalculate the value of the root functions as
                   they cannot be requested from the solver. *)
            discrete s t' (roots_out, roots_out')
          end
      | Cvode.Continue -> continuous s (advtime t')
      | Cvode.StopTimeReached -> finish s t'

  and discrete s t (roots_out, roots_out') =
    (* DISCRETE CALL *)
    (* INV: exists i. roots_in[i] = true *)
    if abs_float (t -. !last_discrete_t) < !zeroc_time_tol
    then incr discrete_tally;

    last_discrete_t := t;

    if !zeroc_limit > 0 && !zeroc_limit <= !discrete_tally then
      (if !zeroc_ignore_at_limit
       then continuous s (advtime t)    (*  ignore  *)
       else raise TooManyZeroCrossings) (* complain *)
    else begin
      if !log then begin
        Cvode.print_time ("R : ", " ") t;
        if !show_root_names
        then
          (Roots.appi (fun i r -> if r then printf "\t%s" roots.(i)) roots_in;
           print_newline ())
        else Roots.print roots_in
      end;
      if lf false roots_in cstates cder roots_out' then begin
        if !log then begin
          print_string "D : ";
          Carray.print_with_time t cstates
        end;
        Carray.clamp !zeroc_deadzone roots_out';
        calculate_roots_in roots_out roots_out';

        if (allow_delta && Roots.exists roots_in)
        then discrete s t (roots_out', roots_out) (* NB: order swapped *)
        else begin
          Cvode.reinit s t cstates;
          Roots.reset roots_in;
          continuous s (advtime t)
        end
      end
      else finish s t
    end

  and finish s t =
    if !log then print_string "--+\n";
    Cvode.free s

  in
  init ()

let run_delta = run true
let run_synchronous = run false

let set_solver l i () =
  lmm := l;
  iter := i

let sprintf = Printf.sprintf

let args n_eq =
  let neq = n_eq - 1 in 
  [
    ("-functional",
     Arg.Unit (set_solver Cvode.Adams Cvode.Functional),
     "(Adams, Functional)");

    ("-dense",
     Arg.Unit (set_solver Cvode.BDF (Cvode.Newton Cvode.Dense)),
     "(BDF, Dense)");

    ("-band",
     Arg.Unit (set_solver Cvode.BDF (Cvode.Newton
        (Cvode.Band { Cvode.mupper = neq; Cvode.mlower = neq }))),
     sprintf "(BDF, Band(%d, %d))" neq neq);

    ("-diag",
     Arg.Unit (set_solver Cvode.BDF (Cvode.Newton Cvode.Diag)),
     "(BDF, Diag)");

    ("-spgmr",
     Arg.Unit (set_solver Cvode.BDF
        (Cvode.Newton (Cvode.Spgmr { Cvode.pretype = Cvode.PrecNone;
                                     Cvode.maxl = 0 }))),
     "(BDF, SPGMR(Both))");

    ("-spbcg",
     Arg.Unit (set_solver Cvode.BDF
        (Cvode.Newton (Cvode.Spbcg { Cvode.pretype = Cvode.PrecNone;
                                     Cvode.maxl = 0 }))),
     "(BDF, SPBCG(Both))");

    ("-sptfqmr",
     Arg.Unit (set_solver Cvode.BDF
        (Cvode.Newton (Cvode.Sptfqmr { Cvode.pretype = Cvode.PrecNone;
                                       Cvode.maxl = 0 }))),
     "(BDF, SPTFQMR(Both))");

    ("-banded-spgmr",
     Arg.Unit (set_solver Cvode.BDF
        (Cvode.Newton (Cvode.BandedSpgmr ({ Cvode.pretype = Cvode.PrecNone;
                                            Cvode.maxl = 0 },
                                          { Cvode.mupper = neq;
                                            Cvode.mlower = neq})))),
     sprintf "(BDF, SPGMR(Both, %d, %d))" neq neq);

    ("-banded-spbcg",
     Arg.Unit (set_solver Cvode.BDF
         (Cvode.Newton (Cvode.BandedSpbcg ({ Cvode.pretype = Cvode.PrecNone;
                                             Cvode.maxl = 0 },
                                           { Cvode.mupper = neq;
                                             Cvode.mlower = neq})))),
     sprintf "(BDF, SPBCG(Both, %d, %d))" neq neq);

    ("-banded-sptfqmr",
     Arg.Unit (set_solver Cvode.BDF
         (Cvode.Newton (Cvode.BandedSptfqmr ({ Cvode.pretype = Cvode.PrecNone;
                                               Cvode.maxl = 0 },
                                             { Cvode.mupper = neq;
                                               Cvode.mlower = neq})))),
     sprintf "(BDF, SPTFQMR(Both, %d, %d))" neq neq);

    ("-simt",
     Arg.Float (fun m -> max_sim_time := Some m),
     "simulation time");

    ("-minstep",
     Arg.Float (fun m -> min_step_size := Some m),
     "minimum step size");

    ("-maxstep", Arg.Set_float max_step_size, "maximum step size");

    ("-onestep",
     Arg.Unit (fun () -> step := Cvode.one_step),
     "Solve the continuous dynamics with CV_ONE_STEP (default: CV_NORMAL).");

    ("-ignore_too_many_zeroc",
     Arg.Set zeroc_ignore_at_limit,
     "Ignore zero-crossings beyond the limit (default: stop simulation).");

    ("-zeroc-limit",
     Arg.Set_int zeroc_limit,
     "Zero-crossing limit (0 = no limit).");

    ("-zeroc-thres",
     Arg.Set_float zeroc_deadzone,
     "Zero-crossing threshold for clamping expression values to zero.");

     ("-reltol",
      Arg.Float (fun t -> rel_tol := Some t),
      "Set relative tolerance (only effective if -abstol is also given).");

     ("-abstol",
      Arg.Float (fun t -> abs_tol := Some t),
      "Set absolute tolerance (only effective if -reltol is also given).");

    ("-precisetime",
     Arg.Set Cvode.extra_time_precision,
     "Plot time values with higher precision.");

    ("-rootnames",
     Arg.Set show_root_names,
     "Show the names of active zero-crossings.");

    ("-addepsilons",
     Arg.Set_int epsilons_to_add,
     "Shift each initial state by the given number of 'epsilons'.");

    ("-l",
     Arg.Set log,
     "Log state variables and zero-crossings to stdout.");
]

let float_with_delta_of_string s =
  let tally = ref 0 in

  let count_deltas c =
    if c == '+' then incr tally
    else if c == '-' then decr tally
  in

  let f = Scanf.sscanf s "%e%s" (fun f s -> (String.iter count_deltas s; f))
  in
  add_epsilons !tally f

let set_float_delta fr =
  Arg.String (fun s -> fr := float_with_delta_of_string s)

