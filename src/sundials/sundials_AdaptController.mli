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

(** Adaptivity controller types and operations.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria) *)

(** Accuracy-based adaptivity controllers for use with various integrators.

    @adaptcontroller SUNAdaptController *)
type t

(** Defines possible controller types.

    @adaptcontroller SUNAdaptController_Type *)
type control_type =
  | NoControl
    (** Performs no control. {cconst SUN_ADAPT_CONTROLLER_NONE} *)
  | SingleRate
    (** Controls a single-rate step size. {cconst SUN_ADAPT_CONTROLLER_H} *)

(** Returns the control type identifier of the controller.

    @adaptcontroller SUNAdaptController_GetType *)
val get_type       : t -> control_type

(** Estimates a single-rate step size.
    The arguments are
    - [c], an adapt controller,
    - [h], the step size from the previous step attempt,
    - [p], the current order of accuracy for the time integration method, and
    - [dsm], the local temporal estimate from the previous step attempt.
    The function returns the estimated step size.

    @adaptcontroller SUNAdaptController_EstimateStep *)
val estimate_step  : t -> h:float -> p:int -> dsm:float -> float

(** Resets a controller to its initial state.

    @adaptcontroller SUNAdaptController_Reset *)
val reset          : t -> unit

(** Sets the controller parameters to their default values.

    @adaptcontroller SUNAdaptController_SetDefaults *)
val set_defaults   : t -> unit

(** Writes all controller parameters to a file (defaults to stdout).

    @adaptcontroller SUNAdaptController_Write *)
val write          : ?logfile:Sundials_impl.Logfile.t -> t -> unit

(** Sets an error bias factor for scaling the local error factors.
    If the given error bias factor is <= 0, the controller uses a
    default value.

    @adaptcontroller SUNAdaptController_SetErrorBias *)
val set_error_bias : t -> float -> unit

(** Notifies a {{!control_type}SingleRate} controller that a
    successful time step was taken with the step size [h] and local error
    factor [dsm].

    @adaptcontroller SUNAdaptController_SetErrorBias *)
val update_h       : t -> h:float -> dsm:float -> unit

(** Returns, respectively, the number of floats and integer values stored
    in the controller.

    @adaptcontroller SUNAdaptController_Space *)
val space          : t -> int * int

(** The Soderlind adaptivity controller implements a general structure for
    temporal control proposed by G. Soderlind.

    @adaptcontroller <SUNAdaptController_links.html#the-sunadaptcontroller-soderlind-module> The SUNAdaptController_Soderlind Module
*)
module Soderlind :
  sig (* {{{ *)

    (** Specify the parameters for a Soderlind adaptivity controller.

        @adaptcontroller SUNAdaptController
        @adaptcontroller SUNAdaptController_Soderlind
        @adaptcontroller SUNAdaptController_SetParams_Soderlind
        @adaptcontroller SUNAdaptController_PID
        @adaptcontroller SUNAdaptController_SetParams_PID
        @adaptcontroller SUNAdaptController_PI
        @adaptcontroller SUNAdaptController_SetParams_PI
        @adaptcontroller SUNAdaptController_I
        @adaptcontroller SUNAdaptController_SetParams_I
        @adaptcontroller SUNAdaptController_ExpGus
        @adaptcontroller SUNAdaptController_SetParams_ExpGus
        @adaptcontroller SUNAdaptController_ImpGus
        @adaptcontroller SUNAdaptController_SetParams_ImpGus *)
    type params =
      | Soderlind of { k1 : float; k2 : float; k3 : float; k4 : float; k5 : float; }
        (** Specify all five Soderlind parameters. *)
      | Defaults
        (** Default Soderlind parameters. *)
      | PID       of { k1 : float; k2 : float; k3 : float; }
        (** Replicate a PID controller ($k_4 = k_5 = 0$). *)
      | PID_defaults
        (** Replicate a PID controller with default parameters. *)
      | PI        of { k1 : float; k2 : float; }
        (** Replicate a PI controller ($k_3 = k_4 = k_5 = 0$). *)
      | PI_defaults
        (** Replicate a PI controller with default parameters. *)
      | I         of { k1 : float; }
        (** Replicate an I controller ($k_2 = k_3 = k_4 = k_5 = 0$). *)
      | I_defaults
        (** Replicate an I controller with default parameters. *)
      | ExpGus    of { k1_hat : float; k2_hat : float; }
        (** Replicate Gustafsson's explicit controller
            ({% $k_1 = \hat{k_1} + \hat{k_2}$ %}, {% $k_2 = -\hat{k_2} %}$,
             {% $k_3 = k_4 = k_5 = 0$ %}). *)
      | ExpGus_defaults
        (** Replicate Gustafsson's explicit controller with default
            parameters. *)
      | ImpGus    of { k1_hat : float; k2_hat : float; }
        (** Replicate Gustafsson's implicit controller
            ({% $k_1 = \hat{k_1} + \hat{k_2}$ %}, {% $k_2 = -\hat{k_2}$ %},
             {% $k_4 = 1$ %}, {% $k_3 = k_5 = 0 $ %}). *)
      | ImpGus_defaults
        (** Replicate Gustafsson's implicit controller with default
            parameters. *)

    (** Create a Soderlind adaptivity controller.
        Uses the default parameters for a Soderline controller, unless
        [~params] is given, in which the appropriate type of controller
        is created.

        @adaptcontroller SUNAdaptController_Soderlind
        @adaptcontroller SUNAdaptController_SetParams_Soderlind
        @adaptcontroller SUNAdaptController_PID
        @adaptcontroller SUNAdaptController_SetParams_PID
        @adaptcontroller SUNAdaptController_PI
        @adaptcontroller SUNAdaptController_SetParams_PI
        @adaptcontroller SUNAdaptController_I
        @adaptcontroller SUNAdaptController_SetParams_I
        @adaptcontroller SUNAdaptController_ExpGus
        @adaptcontroller SUNAdaptController_SetParams_ExpGus
        @adaptcontroller SUNAdaptController_ImpGus
        @adaptcontroller SUNAdaptController_SetParams_ImpGus
        @since 6.7.0 *)
    val make : ?context:Sundials_impl.Context.t -> params -> t

  end (* }}} *)

(** The Gustafsson adaptivity controller implements a combination of two
    adaptivity controllers proposed by K. Gustafsson.

    @adaptcontroller <SUNAdaptController_links.html#the-sunadaptcontroller-imexgus-module> The SUNAdaptController_ImExGus Module
*)
module ImExGus :
  sig (* {{{ *)

    (** Specify the parameters for a Gustafsson adaptivity controller.

        @adaptcontroller SUNAdaptController
        @adaptcontroller SUNAdaptController_ImExGus
        @adaptcontroller SUNAdaptController_SetParams_ImExGus *)
    type params = {
        k1e : float; (** explicit controller parameter 1 *)
        k2e : float; (** explicit controller parameter 2 *)
        k1i : float; (** implicit controller parameter 1 *)
        k2i : float; (** implicit controller parameter 2 *)
      }

    (** Create a Gustafsson implicit-explicit adaptivity controller.
        Uses the default parameters unless [~params] is given.

        @adaptcontroller SUNAdaptController
        @adaptcontroller SUNAdaptController_ImExGus
        @adaptcontroller SUNAdaptController_SetParams_ImExGus
        @since 6.7.0 *)
    val make : ?context:Sundials_impl.Context.t -> params option -> t

  end (* }}} *)

(** Custom implementations of adaptivity controllers. *)
module Custom :
  sig (* {{{ *)

    (** The callback functions required by a custom adaptivity controller.
        The ['s] argument is available for representing an internal state. *)
    type 's adapt_controller_ops = { (* {{{ *)
      get_type       : 's -> control_type;
        (** Return the type of controller. Exceptions are ignored. *)
      estimate_step  : 's -> h:float -> p:int -> dsm:float -> float;
        (** Return the estimated step size given the previous step size [h],
            the current order of accuracy [p], and the previous local temporal
            esimate [dsm]. *)
      reset          : 's -> unit;
        (** Reset to initial state. *)
      set_defaults   : 's -> unit;
        (** Set parameters to their default values. *)
      write          : 's -> Sundials_impl.Logfile.t -> unit;
        (** Write all control parameters to the given file. *)
      set_error_bias : 's -> float -> unit;
        (** Accept an error bias factor for scaling local error factors. *)
      update_h       : 's -> h:float -> dsm:float -> unit;
        (** A successful time step was taken with stepsize [h] and local
            error factor [dsm]. *)
      space          : 's -> int * int;
        (** Return the number of, respectively, floating-point and integer
            values used by the controller. *)
    } (* }}} *)

    (** Create a custom adaptivity controller from a set of callback functions
        and their internal state.

        @adaptcontroller SUNAdaptController
        @since 6.7.0 *)
    val make : ?context:Sundials_impl.Context.t -> 's adapt_controller_ops -> 's -> t

  end (* }}} *)

(** An illegal input was given {cconst SUNADAPTCONTROLLER_ILL_INPUT}.

    @adaptcontroller <SUNAdaptController_links.html#sunadaptcontroller-error-codes> SUNAdaptController Error Codes *)
exception IllInput

(** A user-supplied function failed {cconst SUNADAPTCONTROLLER_USER_FCN_FAIL}.

    @adaptcontroller <SUNAdaptController_links.html#sunadaptcontroller-error-codes> SUNAdaptController Error Codes *)
exception UserFunctionFailure

(** An general error occurred {cconst SUNADAPTCONTROLLER_OPERATION_FAIL}.

    @adaptcontroller <SUNAdaptController_links.html#sunadaptcontroller-error-codes> SUNAdaptController Error Codes *)
exception OperationFailure

