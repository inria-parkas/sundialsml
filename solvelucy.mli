
type solvemode =
| Init          (* set the initial continuous state values *)
| Discrete      (* handle zero-crossings *)
| Continuous    (* solve discrete states *)

type lucyf =
   solvemode
  -> Cvode.Roots.t       (* solvemode = Discrete
                           IN: zero crossings
                          *)
  -> Cvode.val_array     (* solvemode = Init:
                           OUT: initial continuous state values

                           solvemode = Discrete
                           OUT: discrete changes to continuous state values

                           solvemode = Continuous:
                           IN: continuous state values
                          *)
  -> Cvode.der_array     (* solvemode = Continuous
                           OUT: continous derivatives, may be empty
                          *)
  -> Cvode.rootval_array (* solvemode = Continuous
                           OUT: values used for detecting roots

                           solvemode = Discrete
                           OUT: values used for detecting roots
                           NB:  these must be calculated against
                                the updated continuous state values,
                                i.e. not those given to lucyf, but
                                those recalculated within lucyf.
                          *)
  -> bool                (* solvemode = Init or Continuous
                             ignored
                           solvemode = Discrete
                             true: keep running
                             false: quit simulation
                          *)

val sundialify :
  float option ->           (* stop time *)
  lucyf ->                  (* model function *)
  (float -> float) ->       (* advance time *)
  int ->                    (* number of continuous states *)
  int ->                    (* number of zero-crossing functions *)
  unit

