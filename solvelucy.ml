
module Cvode = Cvode_serial;
module Roots = Cvode.Roots


(*
let f zeros conts inputs (* out: *) derivs outputs roots =
  ()
*)

type inputs = int array
let no_inputs = [| |] : int array
(* XXX: Notes:
   - How should we type the discrete inputs?
   - Furthermore, how do we interface changes in discrete inputs with the
     simulation. For example, if the user presses a key / clicks on a block,
     do we signal this as some kind of `external zero crossing'?
 *)
type outputs = int array
let no_outputs = [| |] : int array
(* XXX: Ditto *)

type solvemode =
| Init          (* set the initial continuous state values *)
| Discrete      (* handle zero-crossings *)
| Continuous    (* solve discrete states *)

type lucyf :
   solvemode
  -> Roots.t            (* solvemode = Discrete
                           IN: zero crossings
                         *)
  -> Cvode.val_array    (* solvemode = Init:
                           OUT: initial continuous state values

                           solvemode = Discrete
                           OUT: discrete changes to continuous state values

                           solvemode = Continuous:
                           IN: continuous state values
                         *)
  -> inputs             (* solvemode = Discrete
                           IN: discrete inputs
                         *)
  -> Cvode.der_array    (* solvemode = Continuous
                           OUT: continous derivatives, may be empty
                         *)
  -> outputs            (* solvemode = Discrete
                           OUT: discrete outputs
                         *)
  -> rootval_array      (* solvemode = Continuous
                           OUT: values used for detecting roots

                           solvemode = Discrete
                           OUT: values used for detecting roots
                           NB:  these must be calculated against
                                the updated continuous state values,
                                i.e. not those given to lucyf, but
                                those recalculated within lucyf.
                         *)

(* XXX (related to discussion below):
   - Do we need a special zero-crossing to mark when external inputs have
     occurred?
     Or are they sampled against an internal clock.

   - Likewise, for outputs, should we say when we want them?
     Or rely on internal details.

   - It is the (discrete-)solver's responsibility to recalculate the Roots.t
     input (from the Cvode.val_array).
 *)



(* Compare roots_in to roots_in' and return:
   true -> new zero-crossings have occurred
   false -> no new zero-crossings have occurred
 *)
let roots_policy roots_in roots_in' =
  let rin  = Roots.get roots_in
  and rin' = Roots.get roots_in'
  in
  let check i =
    if i < 0 then return false
    else if rin i <> rin' i then true
    else check (i - 1)
  in
  check ((min (Roots.length roots_in) (Roots.length roots_in')) - 1)

(* calculate ri by comparing ro (before) to ro (after). *)
let calculate_roots ri ro ro' =
  let rin = Roots.set ri
  let f i =
    if i < 0 then ()
    else begin
      rin i (ro.{i} < 0 && ro'.{i} >= 0); (* TODO: check definition *)
      f (i - 1)
    end
  in
  f (n_roots - 1)

(*
 * TODO:
 * Currently sundialify runs indefinitely.
 * Should it return from time to time (with frozen state) to allow the caller to
 * update inputs and outputs, and to draw simulation frames, log data, etc?
 * How should this interface work?

   Related questions:
   - Issue: how to handle stopping the simulation.
            Take a final time as an argument?
            Have a special discrete output (terminate)?

   - Issue: how to handle the display of results?
            extra functions?
            do it within the lucyf?

 *)
let sundialify (lf : lucyf) (advtime : float -> float) (n_cstates, n_roots) inputs outputs =
  let cstates    = Cvode.create n_cstates
  and cder       = Cvode.create n_cstates

  and roots_in   = Roots.create n_roots
  and roots_in'  = Roots.create n_roots

  and roots_out  = Cvode.create n_roots
  and roots_out' = Cvode.create n_roots
  in

  let f t cs ds =
    lf Continuous Roots.empty cs no_inputs ds no_outputs Roots.empty;
  and g t cs rs =
    lf Continuous Roots.empty cs no_inputs Cvode.empty no_outputs rs;
  in

  let calculate_roots_out () = g t' cstates roots_out

  let rec init () =
    Cvode.fill cder 0.0;
    Roots.reset roots_in;

    (* INIT CALL *)
    lf Init Roots.empty cstates no_inputs CVode.empty no_outputs Roots.empty;

    let s = Cvode.init Cvode.Adams Cvode.Functional f (n_roots, g) cstates
    in continuous s (advtime 0.0)

  and continuous s t =
    (* CONTINUOUS CALL(S) *)
    let (t', roots_found) = Cvode.advance s t cstates in
    if roots_found
    then begin
      Cvode.get_roots s roots_in;
      calculate_roots_out ();
      (* NB: we are forced to recalculate the value of the root functions as
             they cannot be requested from the solver. *)
      discrete s t' (roots_out, roots_out')
    end
    else begin
      continuous s (advtime t')
    end

  and discrete s t (roots_out, roots_out') =
    (* DISCRETE CALL *)
    lf Discrete roots_in cstates inputs CVode.empty outputs roots_out';
    calculate_roots_in roots_in' roots_out roots_out';

    if roots_policy roots_in roots_in'
    then discrete s t (roots_out', roots_out) (* NB: order swapped *)
    else begin
      Cvode.reinit s t cstates;
      continuous s (advtime t)
    end
  in
  init ()

(* TODO:
   - Possible problem: advtime gives an increment that is not big enough for the
     solver. So catch exception Cvode.TooClose (inside continuous) and try again?

   - Setup the solver so that it only looks for positive-going roots, per the
     default for lucy-h.
 *)


(*
   XXX Notes to Tim:
   - There is no difference between continuous outputs and continuous states (Moore = Mealy)
   - But there is a difference between discrete outputs and discrete states (Moore /= Mealy)
     and, in fact, last memories conflate the two.
     whereas flows do not.

 *)
