(* [sundials 2.5.0]

   This test shows a problem with Ida.set_stop_time.  The stop time is checked
   to ensure you're not setting a time that has passed, but IDA does this check
   against its current *internal* time, not the last time you asked the value
   at.  So even if your last query to Ida.solve_normal was t, there's no
   guarantee that set_stop_time (t +. 1.) would succeed.

   The only ways to avoid this problem are:

   - set stop time before any calls to solve_normal, solve_one_step, calc_ic_y,
     or calc_ic_ya_yd'; or

   - check where the solver's internal clock is at, using get_current_time

   Also, the value of "current t" in the error printed out when the check fails
   is garbage (too few arguments to vsnprintf() causes invalid access).

 *)

module Ida = Ida_serial;;
module Carray = Ida.Carray;;
let vec = Carray.of_array [| 0. |] and vec' = Carray.of_array [| 1. |]
let session =
  Ida.init Ida.Dense
    (fun t vec vec' res -> res.{0} <- vec'.{0})
    ~t0:0. vec vec'
let _ = Ida.solve_normal session 1. vec vec'
let _ = Ida.set_stop_time session 1.25
