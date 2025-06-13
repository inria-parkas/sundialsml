(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2025 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

exception IllInput
exception UserFunctionFailure
exception OperationFailure

type t

(* synchronized with sundials_ml.h: sundials_adapt_soderlind_control_type_tag *)
type control_type =
  | NoControl
  | SingleRate

external get_type : t -> control_type
  = "sunml_adapt_get_type"

external c_estimate_step : t -> float -> int -> float -> float
  = "sunml_adapt_estimate_step"

let estimate_step s ~h ~p ~dsm = c_estimate_step s h p dsm

external reset : t -> unit
  = "sunml_adapt_reset"

external set_defaults : t -> unit
  = "sunml_adapt_set_defaults"

external c_write : t -> Sundials_impl.Logfile.t -> unit
  = "sunml_adapt_write"

let write ?(logfile=Sundials_impl.Logfile.stdout) s = c_write s logfile

external set_error_bias : t -> float -> unit
  = "sunml_adapt_set_error_bias"

external c_update_h : t -> float -> float -> unit
  = "sunml_adapt_update_h"

let update_h s ~h ~dsm = c_update_h s h dsm

external space : t -> int * int
  = "sunml_adapt_space"

module Soderlind =
  struct (* {{{ *)

    (* synchronized with sundials_ml.h: sundials_adapt_soderlind_params_tag *)
    type params =
      | Soderlind of { k1 : float; k2 : float; k3 : float; k4 : float; k5 : float; }
      | Defaults
      | PID       of { k1 : float; k2 : float; k3 : float; }
      | PID_defaults
      | PI        of { k1 : float; k2 : float; }
      | PI_defaults
      | I         of { k1 : float; }
      | I_defaults
      | ExpGus    of { k1_hat : float; k2_hat : float; }
      | ExpGus_defaults
      | ImpGus    of { k1_hat : float; k2_hat : float; }
      | ImpGus_defaults

    external c_make : Sundials_impl.Context.t -> params -> t
      = "sunml_adapt_soderlind_make"

    let make ?context params =
      let ctx = Sundials_impl.Context.get context in
      c_make ctx params

  end (* }}} *)

module ImExGus =
  struct (* {{{ *)

    type params = {
        k1e : float;
        k2e : float;
        k1i : float;
        k2i : float;
      }

    external c_make : Sundials_impl.Context.t -> params option -> t
      = "sunml_adapt_imexgus_make"

    let make ?context params =
      let ctx = Sundials_impl.Context.get context in
      c_make ctx params

  end (* }}} *)

module Custom =
  struct (* {{{ *)

    type 's adapt_controller_ops = {
      get_type       : 's -> control_type;
      estimate_step  : 's -> h:float -> p:int -> dsm:float -> float;
      reset          : 's -> unit;
      set_defaults   : 's -> unit;
      write          : 's -> Sundials_impl.Logfile.t -> unit;
      set_error_bias : 's -> float -> unit;
      update_h       : 's -> h:float -> dsm:float -> unit;
      space          : 's -> int * int;
    }

    external c_make : Sundials_impl.Context.t -> ('s * 's adapt_controller_ops) -> t
      = "sunml_adapt_custom_make"

    let make ?context ops s =
      let ctx = Sundials_impl.Context.get context in
      c_make ctx (s, ops)

  end (* }}} *)

