(* [sundials 2.5.0]

   This test shows a problem with Cvode.set_stop_time.

   CVODE takes a step in the wrong direction and yields incorrect values when
   stop time was set in the past before the first call to
   Cvode.solve_normal. or Cvode.solve_one_step.  Note that since stop times
   persist across reinits, this means that if you use reinit to set t0 to a
   time past the stop time, then the next call to solve_normal or
   solve_one_step will fail.

   The comment in CVodeSetStopTime() says CVODE is meant to detect and report
   stop times being in the past, but the code that does that check is buggy.
   This line in CVode():

      if ( tstopset && (tout-tn)*(tout-tstop) > 0 ) tout_hin = tstop;

   sets the goal time to tstop when tout,tn are on the same side of tout, i.e.
   when

      tn < tstop < tout,
      tstop < tn < tout,
      tout < tn < tstop, or
      tout < tstop < tn

   where

      tn    = internal time of the session (= t0 in the case of our interest)
      tout  = the time to integrate to (argument passed to CVode())
      tstop = stop time.

   Let's say we're integrating forward, so tn < tout.  If tstop < tout as well,
   then CVODE integrates in the direction of tstop and interpolates tout.  If
   it happens that tstop < tn < tout, that direction is negative.  So CVODE
   integrates from tn down to tstop, tries to interpolate at tout when in fact
   tstop < tn < tout means it's extrapolation, and fails.

   Additionally, CVODE ignores the error return from the interpolator, so it
   reports successful integration to tout and returns an incorrect vector.

 *)

module Cvode = Cvode_serial;;
module Carray = Cvode.Carray;;
let vec = Carray.of_array [| 0. |]
let session =
  Cvode.init Cvode.BDF Cvode.Functional
    (fun t vec vec' -> vec'.{0} <- 1.)
    ~t0:0. vec
let _ =
  Cvode.set_stop_time session 1.;       (* set stop time *)
  (* integrate in the other direction *)
  let (tret, _) = Cvode.solve_normal session (-1.) vec in
  (* we get incorrect values *)
  Printf.printf "tret = %g, vec = [|%g|]\n" tret vec.{0}
