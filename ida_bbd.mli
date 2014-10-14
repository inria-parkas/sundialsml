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
(*               User Documentation for IDA v2.7.0                     *)
(*                Alan C. Hindmarsh and Radu Serban                    *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(** Parallel band-block-diagonal preconditioners for IDA (requires MPI).

    @version VERSION()
    @author Timothy Bourke (Inria)
    @author Jun Inoue (Inria)
    @author Marc Pouzet (LIENS)
    @ida <node5#sss:idabbdpre> Parallel band-block-diagonal preconditioner module *)

(** An alias for sessions based on parallel nvectors. *)
type parallel_session =
      (Nvector_parallel.data, Nvector_parallel.kind) Ida.session

(** An alias for preconditioners based on parallel nvectors. *)
type parallel_preconditioner =
      (Nvector_parallel.data, Nvector_parallel.kind) Ida.Spils.preconditioner

(** The bandwidths for the difference quotient Jacobian operation. *)
type bandwidths = Ida_impl.IdaBbdTypes.bandwidths =
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

(** [gloc t y gd] computes [g(t, y)] into [gd].
    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @idas <node5#sss:idabbdpre> IDABBDLocalFn *)
type local_fn = float
                -> Nvector_parallel.data
                -> Nvector_parallel.data
                -> Nvector_parallel.data
                -> unit

(** [cfn t y] performs all interprocess communication necessary
    for the execution of [local_fn] using the input vector [y].
    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @idas <node5#sss:idabbdpre> IDABBDCommFn *)
type comm_fn = float -> Nvector_parallel.data -> Nvector_parallel.data -> unit

(** Same as {!Ida.Spils.prec_left} but uses the Parallel
    Band-Block-Diagonal preconditioner included in IDA.  Called like
    [prec_left ~dqrely:dqrely bandwidths callbacks], where:

    - [~dqrely] gives the relative increment in components of [y] used in
      the difference quotient approximations
      (defaults to [sqrt unit_roundoff]).
    - [bandwidths] specify the bandwidths to be used in the difference
      quotient Jacobian operation.
    - [callbacks] gives the preconditioning callbacks.  See the
      {!callbacks} type.

    @ida <node5#sss:lin_solv_init> IDASpgmr
    @ida <node5#sss:idabbdpre> IDABBDPrecInit *)
val prec_left : ?dqrely:float
                -> bandwidths
                -> ?comm_fn:comm_fn
                -> local_fn
                -> parallel_preconditioner

(** [reinit s mudq mldq ~dqrely:dqrely] reinitializes the BBD preconditioner
    with upper ([mudq]) and lower ([mldq]) half-bandwidths to be used in the
    difference quotient Jacobian approximation, and an optional relative
    increment in components of [y] (passing [None] uses the default value [sqrt
    unit_roundoff]).

    @ida <node5#sss:idabbdpre> IDABBDPrecReInit *)
val reinit : parallel_session -> ?dqrely:float -> int -> int -> unit

(** {4 Optional output functions} *)

(** Returns the sizes of the real and integer workspaces used by the
    band-block-diagonal preconditioner module.

    @ida <node5#sss:idabbdpre> IDABBDPrecGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space : parallel_session -> int * int

(** Returns the number of calls made to the user-supplied right-hand
    side function due to finite difference banded Jacobian approximation in the
    preconditioner setup function.

    @ida <node5#sss:idabbdpre> IDABBDPrecGetNumGfnEvals *)
val get_num_gfn_evals : parallel_session -> int

