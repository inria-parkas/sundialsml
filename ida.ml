(***********************************************************************)
(*                                                                     *)
(*     OCaml interface to Sundials (serial) CVODE and IDA solvers      *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2013 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

(*
 * NB: The order of variant constructors and record fields is important!
 *     If these types are changed or augmented, the corresponding declarations
 *     in cvode_serial.h (and code in cvode_serial.c) must also be updated.
 *)

include Sundials
include Dls

type root_direction = RootDirs.root_direction

(* Solver exceptions *)
exception IllInput
exception TooClose
exception TooMuchWork
exception TooMuchAccuracy
exception ErrFailure
exception ConvergenceFailure
exception LinearInitFailure
exception LinearSetupFailure
exception LinearSolveFailure
exception ResFuncFailure
exception FirstResFuncFailure
exception RepeatedResFuncErr
exception UnrecoverableResFuncErr
exception RootFuncFailure

(* Initial condition calculator exceptions *)
exception NoRecovery
exception BadEwt

(* get_dky exceptions *)
exception BadK
exception BadT
exception BadDky

let no_roots = (0, (fun _ _ _ _ -> ()))

(* Throw inside the f callback if the derivatives cannot be calculated at
   the given time. *)
exception RecoverableFailure

type integrator_stats = {
    num_steps : int;
    num_res_evals : int;
    num_lin_solv_setups : int;
    num_err_test_fails : int;
    last_order : int;
    current_order : int;
    actual_init_step : float;
    last_step : float;
    current_step : float;
    current_time : float
  }

exception StopTimeReached

exception ZeroDiagonalElement of int

module VarTypes =
  struct
    type t = Carray.t
    type var_type = Algebraic | Differential

    let var_type_of_float = function
      | 0.0 -> Algebraic
      | 1.0 -> Differential
      | f -> raise (Invalid_argument
                      ("invalid component type: " ^ string_of_float f))
    let float_of_var_type = function
      | Algebraic -> 0.0
      | Differential -> 1.0

    let create = Carray.create
    let init n x = Carray.init n (float_of_var_type x)
    let of_array a =
      let ret = create (Array.length a) in
      for i = 0 to Array.length a - 1 do
        ret.{i} <- float_of_var_type a.(i)
      done;
      ret
    let length = Carray.length

    let get a i = var_type_of_float a.{i}
    let set a i x = a.{i} <- float_of_var_type x
    let set_algebraic a i = set a i Algebraic
    and set_differential a i = set a i Differential
    let fill a t =
      let x = float_of_var_type t in
      Carray.fill a x

    let blit a b = Carray.blit a b
  end
module Id = VarTypes

let _ =
  List.iter (fun (nm, ex) -> Callback.register_exception nm ex)
  [
    ("ida_RecoverableFailure",      RecoverableFailure);

    ("ida_StopTimeReached",         StopTimeReached);
    ("ida_IllInput",                IllInput);
    ("ida_TooClose",                TooClose);
    ("ida_TooMuchWork",             TooMuchWork);
    ("ida_TooMuchAccuracy",         TooMuchAccuracy);
    ("ida_ErrFailure",              ErrFailure);
    ("ida_ConvergenceFailure",      ConvergenceFailure);
    ("ida_LinearInitFailure",       LinearInitFailure);
    ("ida_LinearSetupFailure",      LinearSetupFailure);
    ("ida_LinearSolveFailure",      LinearSolveFailure);
    ("ida_ResFuncFailure",          ResFuncFailure);
    ("ida_FirstResFuncFailure",     FirstResFuncFailure);
    ("ida_RepeatedResFuncErr",      RepeatedResFuncErr);
    ("ida_NoRecovery",              NoRecovery);
    ("ida_BadEwt",                  BadEwt);
    ("ida_UnrecoverableResFuncErr", UnrecoverableResFuncErr);
    ("ida_RootFuncFailure",         RootFuncFailure);

    ("ida_BadK",                    BadK);
    ("ida_BadT",                    BadT);
    ("ida_BadDky",                  BadDky);

    ("ida_ZeroDiagonalElement",     ZeroDiagonalElement 0);
  ]
