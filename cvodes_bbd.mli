(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

(***********************************************************************)
(* Much of the comment text is taken directly from:                    *)
(*                                                                     *)
(*               User Documentation for CVODE v2.6.0                   *)
(*                Alan C. Hindmarsh and Radu Serban                    *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(** Parallel band-block-diagonal preconditioners for CVODES (requires MPI).

    @version VERSION()
    @author Timothy Bourke (Inria)
    @author Jun Inoue (Inria)
    @author Marc Pouzet (LIENS)
    @cvodes <node7#SECTION00742000000000000000> Using the band-block-diagonal preconditioner CVBBDPRE
 *)

(** An alias for sessions based on parallel nvectors. *)
type parallel_bsession =
      (Nvector_parallel.data, Nvector_parallel.kind) Cvodes.Adjoint.bsession

(** An alias for preconditioners based on parallel nvectors. *)
type parallel_preconditioner = (Nvector_parallel.data, Nvector_parallel.kind)
                                            Cvodes.Adjoint.Spils.preconditioner

(** The bandwidths for the difference quotient Jacobian operation. *)
type bandwidths = Cvode_bbd.bandwidths =
  {
    mudq    : int; (** Upper half-bandwidth for the difference
                       quotient Jacobian approximation. *)
    mldq    : int; (** Lower half-bandwidth for the difference
                       quotient Jacobian approximation. *)
    mukeep  : int; (** Upper half-bandwidth for the retained banded
                       approximate Jacobian block. *)
    mlkeep  : int; (** Lower half-bandwidth for the retained banded
                       approximate Jacobian block. *)
  }

(** [gloc t y yb gb] computes [g(t, y)] into [gb] from the value of the
    independent variable [t], the current forward solution vector [y], and
    the current value of the backward dependent variable vector.
    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @cvodes <node7#SECTION00742200000000000000> CVBBDLocalFnB *)
type local_fn = float
                -> Nvector_parallel.data
                -> Nvector_parallel.data
                -> Nvector_parallel.data
                -> unit

(** [cfn t y yb] performs all interprocess communication necessary for the
    execution of [local_fn] using the forward solution vector [y] and the
    backward dependent variable vector [yb].
    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @cvodes <node7#SECTION00742200000000000000> CVBBDCommFnB *)
type comm_fn = float -> Nvector_parallel.data -> Nvector_parallel.data -> unit

(** Same as {!Cvodes.Adjoint.Spils.prec_left} but uses the Parallel
    Band-Block-Diagonal preconditioner included in CVODES.  Called like
    [prec_left ~dqrely:dqrely bandwidths callbacks], where:

    - [~dqrely] gives the relative increment in components of [y] used in
      the difference quotient approximations
      (defaults to [sqrt unit_roundoff]).
    - [bandwidths] specify the bandwidths to be used in the difference
      quotient Jacobian operation.
    - [callbacks] gives the preconditioning callbacks.  See the
      {!callbacks} type.

    @cvodes <node7#SECTION00742100000000000000> CVBBDPrecInitB *)
val prec_left : ?dqrely:float
                -> bandwidths
                -> ?comm_fn:comm_fn
                -> local_fn
                -> parallel_preconditioner

(** Same as {!prec_left} but preconditions from the right.

    @cvodes <node7#SECTION00742100000000000000> CVBBDPrecInitB *)
val prec_right : ?dqrely:float
                 -> bandwidths
                 -> ?comm_fn:comm_fn
                 -> local_fn
                 -> parallel_preconditioner

(** Same as {!prec_left} but preconditions from both sides.

    @cvodes <node7#SECTION00742100000000000000> CVBBDPrecInitB *)
val prec_both : ?dqrely:float
                -> bandwidths
                -> ?comm_fn:comm_fn
                -> local_fn
                -> parallel_preconditioner

(** [reinit s mudq mldq ~dqrely:dqrely] reinitializes the BBD
    preconditioner with upper ([mudq]) and lower ([mldq])
    half-bandwidths to be used in the difference quotient Jacobian
    approximation, and an optional relative increment [dqrely] in
    components of [y].  [dqrely] defaults to
    [sqrt Sundials.unit_roundoff].

    @cvodes <node7#SECTION00742000000000000000> CVBBDPrecReInitB *)
val reinit : parallel_bsession -> ?dqrely:float -> int -> int -> unit

(** {4 Optional output functions} *)

(** Returns the sizes of the real and integer workspaces used by the
    band-block-diagonal preconditioner module.

    @cvode <node7#SECTION00742000000000000000> CVBBDPrecGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space : parallel_bsession -> int * int

(** Returns the number of calls made to the user-supplied right-hand
    side function due to finite difference banded Jacobian approximation in the
    preconditioner setup function.

    @cvode <node7#SECTION00742000000000000000> CVBBDPrecGetNumGfnEvals *)
val get_num_gfn_evals : parallel_bsession -> int

