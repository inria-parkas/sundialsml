
module Cvode = Cvode_serial
module Roots = Cvode.Roots
module Carray = Cvode.Carray

type solvemode =
| Init          (* set the initial continuous state values *)
| Discrete      (* handle zero-crossings *)
| Continuous    (* solve discrete states *)

type lucyf =
   solvemode
  -> Roots.t
  -> Cvode.val_array
  -> Cvode.der_array
  -> Cvode.rootval_array
  -> bool

(* Compare roots_in to roots_in' and return:
   true -> new zero-crossings have occurred
   false -> no new zero-crossings have occurred
 *)
let roots_policy roots_in roots_in' =
  let rin  = Roots.get roots_in
  and rin' = Roots.get roots_in'
  in
  let rec check i =
    if i < 0 then false
    else if rin i <> rin' i then true
    else check (i - 1)
  in
  check ((min (Roots.length roots_in) (Roots.length roots_in')) - 1)

(* calculate ri by comparing ro (before) to ro (after). *)
let calculate_roots ri ro ro' =
  let rin = Roots.set ri in
  let rec f i =
    if i < 0 then ()
    else begin
      rin i (ro.{i} < 0.0 && ro'.{i} >= 0.0); (* TODO: check definition *)
      f (i - 1)
    end
  in
  f ((Carray.length ro) - 1)

let sundialify tmax (lf : lucyf) (advtime : float -> float) (n_cstates, n_roots) =
  let cstates    = Carray.create n_cstates
  and cder       = Carray.create n_cstates

  and roots_in   = Roots.create n_roots
  and roots_in'  = Roots.create n_roots

  and roots_out  = Carray.create n_roots
  and roots_out' = Carray.create n_roots
  in

  let f t cs ds =
    ignore (lf Continuous Roots.empty cs ds Carray.empty);
  and g t cs rs =
    ignore (lf Continuous Roots.empty cs Carray.empty rs);
  in

  let calculate_roots_out t = g t cstates roots_out in

  let rec init () =
    Carray.fill cder 0.0;
    Roots.reset roots_in;

    (* INIT CALL *)
    ignore (lf Init Roots.empty cstates Carray.empty Carray.empty);

    let s = Cvode.init Cvode.Adams Cvode.Functional f (n_roots, g) cstates in
    Cvode.set_all_root_directions s Cvode.Increasing;
    match tmax with None -> () | Some t -> Cvode.set_stop_time s t;
    continuous s (advtime 0.0)

  and continuous s t =
    (* CONTINUOUS CALL(S) *)
    let (t', roots_found) = Cvode.advance s t cstates
    in
      if roots_found
      then begin
        Cvode.get_roots s roots_in;
        calculate_roots_out t';
        (* NB: we are forced to recalculate the value of the root functions as
               they cannot be requested from the solver. *)
        discrete s t' (roots_out, roots_out')
      end
      else begin
        continuous s (advtime t')
      end

  and discrete s t (roots_out, roots_out') =
    (* DISCRETE CALL *)
    if lf Discrete roots_in cstates Carray.empty roots_out' then begin
      calculate_roots roots_in' roots_out roots_out';

      if roots_policy roots_in roots_in'
      then discrete s t (roots_out', roots_out) (* NB: order swapped *)
      else begin
        Cvode.reinit s t cstates;
        continuous s (advtime t)
      end
    end
    else finish s t

  and finish s t =
    Cvode.free s

  in
  init ()

(* TODO:
   - Think harder about the interface between simulation code and the external
     world?
     
     For instance, if a discrete node calls a function to get the mouse
     position, and there are multiple Discrete iterations at an instant, is it
     important to hold this value constant (which may be expensive, in terms of
     memory, and complicated), or just let it make multiple calls?

     Worse, what about a destructive input, like reading bytes from a file; it
     should probably only read one value regardless of the number of
     zero-crossing iterations.

     Even more so for outputs and state changes.

   - Do we need a special zero-crossing to mark when external inputs have
     occurred?
     Or are they sampled against an internal clock; i.e. does the discrete
     program have to specify a sampling rate.
 *)


(*
   XXX Notes to Tim:
   - There is no difference between continuous outputs and continuous states (Moore = Mealy)
   - But there is a difference between discrete outputs and discrete states (Moore /= Mealy)
     and, in fact, last memories conflate the two.
     whereas flows do not.

   - It should be a piece of cake to reimplement the Argos in Simulink model
     using hybrid lucid synchrone.
 *)

