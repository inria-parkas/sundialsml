
type solvemode =
| Init          (* set the initial continuous state values *)
| Discrete      (* handle zero-crossings *)
| Continuous    (* solve discrete states *)

type lucyf =
   solvemode
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

(* In brief, by mode:

   mode == Init
        calculate: y

   mode == Discrete
        using: rin
        calculate: y, rout
        return: true (to continue) or false (to terminate)

   mode == Continuous
        using: y
        calculate: der, rout
 *)

val sundialify :
  float option ->           (* stop time *)
  lucyf ->                  (* model function *)
  (float -> float) ->       (* advance time *)
  int ->                    (* number of continuous states *)
  int ->                    (* number of zero-crossing functions *)
  unit

