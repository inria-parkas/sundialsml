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

(** Parallel band-block-diagonal preconditioners for CVODE (requires MPI).

    @version VERSION()
    @author Timothy Bourke (Inria)
    @author Jun Inoue (Inria)
    @author Marc Pouzet (LIENS)
    @cvode <node5#sss:cvbbdpre> Parallel band-block-diagonal preconditioner module *)

(** An alias for sessions based on parallel nvectors. *)
type parallel_session =
      (Nvector_parallel.data, Nvector_parallel.kind) Cvode.session

(** An alias for preconditioners based on parallel nvectors. *)
type parallel_preconditioner =
      (Nvector_parallel.data, Nvector_parallel.kind) Cvode.Spils.preconditioner

(** The bandwidths for the difference quotient Jacobian operation. *)
type bandwidths = Cvode_impl.CvodeBbdTypes.bandwidths =
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

(* TODO: *)
(** [gloc t y gd] computes [g(t, y)] into [gd].

    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @cvodes <node5#sss:cvbbdpre> CVBBDLocalFn *)
type local_fn = float -> Nvector_parallel.data -> Nvector_parallel.data -> unit

(** Functions that perform the interprocess communication necessary
    for the execution of {!local_fn}.
    In the call [cfn t y], [t] is the independent variable (time) and [y] is
    the input vector.

    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @cvodes <node5#sss:cvbbdpre> CVBBDCommFn *)
type comm_fn = float -> Nvector_parallel.data -> unit

(** Left preconditioning using the Parallel Band-Block-Diagonal module.
    In the call [prec_left ~dqrely:dqrely bandwidths callbacks], the
    difference quotient Jacobian operation is controlled by [dqrely], which
    gives the relative increment in components of [y], and [bandwidths].

    @cvode <node5#sss:cvbbdpre> CVBBDPrecInit *)
val prec_left : ?dqrely:float
                -> bandwidths
                -> ?comm_fn:comm_fn
                -> local_fn
                -> parallel_preconditioner

(** Same as {!prec_left} but preconditions from the right.

    @cvode <node5#sss:cvbbdpre> CVBBDPrecInit *)
val prec_right : ?dqrely:float
                 -> bandwidths
                 -> ?comm_fn:comm_fn
                 -> local_fn
                 -> parallel_preconditioner

(** Same as {!prec_left} but preconditions from both
    sides.

    @cvode <node5#sss:cvbbdpre> CVBBDPrecInit *)
val prec_both : ?dqrely:float
                -> bandwidths
                 -> ?comm_fn:comm_fn
                 -> local_fn
                -> parallel_preconditioner

(** [reinit s mudq mldq ~dqrely:dqrely] reinitializes the BBD
    preconditioner with upper ([mudq]) and lower ([mldq])
    half-bandwidths to be used in the difference quotient Jacobian
    approximation, and an optional relative increment [dqrely] in
    components of [y].  [dqrely] defaults to [sqrt unit_roundoff].

    @cvode <node5#sss:cvbbdpre> CVBBDPrecReInit *)
val reinit : parallel_session -> ?dqrely:float -> int -> int -> unit

(** {4 Optional output functions} *)

(** Returns the sizes of the real and integer workspaces used by the
    band-block-diagonal preconditioner module.

    @cvode <node5#sss:cvbbdpre> CVBBDPrecGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space : parallel_session -> int * int

(** Returns the number of calls made to the user-supplied right-hand
    side function due to finite difference banded Jacobian approximation in the
    preconditioner setup function.

    @cvode <node5#sss:cvbbdpre> CVBBDPrecGetNumGfnEvals *)
val get_num_gfn_evals : parallel_session -> int

