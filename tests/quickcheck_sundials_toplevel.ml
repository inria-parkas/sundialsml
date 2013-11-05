(* This module should be loaded when using the sundials part of quickcheck from
   the top level.  It sets a flag so that certain problems can be avoided.  See
   the comment where inhibit_quickcheck_main is defined.  *)
Quickcheck_sundials.inhibit_quickcheck_main := true
