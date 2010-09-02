
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

let sundialify tmax (lf : lucyf) advtime n_cstates n_roots =
  let cstates    = Carray.create n_cstates
  and cder       = Carray.create n_cstates

  and roots_in   = Roots.create n_roots
  and no_roots_in = Roots.create n_roots

  and roots_out  = Carray.create n_roots
  and roots_out' = Carray.create n_roots
  in

  let f t cs ds =
    ignore (lf false roots_in cs ds roots_out)
  and g t cs rs =
    ignore (lf false no_roots_in cs cder rs)
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

    let s = Cvode.init Cvode.Adams Cvode.Functional f (n_roots, g) cstates in
    Cvode.set_all_root_directions s Cvode.Increasing;
    match tmax with None -> () | Some t -> Cvode.set_stop_time s t;
    Roots.reset roots_in;
    continuous s (advtime 0.0)

  and continuous s t =
    (* CONTINUOUS CALL(S) *)
    (* INV: forall i. roots_in[i] = false *)
    let (t', result) = Cvode.normal s t cstates
    in
      print_string "C: "; (* XXX *)
      Carray.print_with_time t' cstates; (* TODO: how to handle display in general *)
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
    print_string "R: "; Roots.print roots_in; (* TODO: how to handle display in general *)
    if lf false roots_in cstates cder roots_out' then begin
      print_string "D: "; (* XXX *)
      Carray.print_with_time t cstates; (* TODO: how to handle display in general *)
      calculate_roots_in roots_out roots_out';

      if Roots.exists roots_in
      then discrete s t (roots_out', roots_out) (* NB: order swapped *)
      else begin
        Cvode.reinit s t cstates;
        Roots.reset roots_in;
        continuous s (advtime t)
      end
    end
    else finish s t

  and finish s t =
    Cvode.free s

  in
  init ()

