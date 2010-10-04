(*
 * Timothy Bourke (INRIA) & Marc Pouzet (ENS), August 2009
 *
 * Prototype solver for hybrid-lucy programs compiled to a single function.
 * Implements both the 'delta-step' semantics (run_delta), i.e. instantaneous
 * loop on new zero-crossings in the discrete step, and the 'synchronous'
 * semantics (run_synchronous), i.e. a single execution at each discrete step.
 *)

type lucyf =
   bool                         (* true: init, false: continuous/discrete *)
  -> Cvode_serial.Roots.t       (* solvemode = Discrete
                                     IN: zero crossings
                                *)
  -> Cvode_serial.val_array     (* solvemode = Init:
                                     OUT: initial continuous state values

                                     solvemode = Discrete
                                     OUT: discrete changes to continuous state values

                                     solvemode = Continuous:
                                     IN: continuous state values
                                *)
  -> Cvode_serial.der_array     (* solvemode = Continuous
                                     OUT: continous derivatives, may be empty
                                *)
  -> Cvode_serial.rootval_array (* solvemode = Continuous
                                     OUT: values used for detecting roots

                                   solvemode = Discrete
                                     OUT: values used for detecting roots
                                     NB:  these must be calculated against
                                          the updated continuous state values,
                                          i.e. not those given to lucyf, but
                                          those recalculated within lucyf.
                                *)
  -> bool                       (* solvemode = Init or Continuous
                                     ignored
                                   solvemode = Discrete
                                     true: keep running
                                     false: quit simulation
                                *)

(* The mode is determined by:
 *   init = true --> mode == Init
 *   init = false && forall i. roots[i] == false --> mode == Continuous
 *   init = false && exists i. roots[i] == true  --> mode == Discrete
 *
 * In brief, by mode:
 *
 * mode == Init
 *      calculate: y
 *
 * mode == Discrete
 *      using: rin
 *      calculate: y, rout
 *      return: true (to continue) or false (to terminate)
 *
 * mode == Continuous
 *      using: y
 *      calculate: der, rout
 *)

val lmm : Cvode_serial.lmm ref
val iter : Cvode_serial.iter ref
val enable_logging : unit -> unit
val enable_zeroc_logging : unit -> unit

val args : int -> (Arg.key * Arg.spec * Arg.doc) list
val set_float_delta : float ref -> Arg.spec

val run :
  bool ->                   (* allow multiple discrete delta-steps *)
  float option ->           (* stop time *)
  lucyf ->                  (* model function *)
  (float -> float) ->       (* advance time *)
  int ->                    (* number of continuous states *)
  int ->                    (* number of zero-crossing functions *)
  unit

val run_delta :
  float option ->           (* stop time *)
  lucyf ->                  (* model function *)
  (float -> float) ->       (* advance time *)
  int ->                    (* number of continuous states *)
  int ->                    (* number of zero-crossing functions *)
  unit

val run_synchronous :
  float option ->           (* stop time *)
  lucyf ->                  (* model function *)
  (float -> float) ->       (* advance time *)
  int ->                    (* number of continuous states *)
  int ->                    (* number of zero-crossing functions *)
  unit

