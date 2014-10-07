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

type data = Nvector_parallel.data
type kind = Nvector_parallel.kind
type parallel_session = (data, kind) Ida.session
type parallel_preconditioner = (data, kind) Ida.Spils.preconditioner

type bandwidths = Ida_impl.IdaBbdTypes.bandwidths =
  {
    mudq    : int; (** Upper half-bandwidth to be used in the difference
                       quotient Jacobian approximation. *)
    mldq    : int; (** Lower half-bandwidth to be used in the difference
                       quotient Jacobian approximation. *)
    mukeep  : int; (** Upper half-bandwidth of the retained banded approximate
                       Jacobian block. *)
    mlkeep  : int; (** Lower half-bandwidth of the retained banded approximate
                       Jacobian block. *)
  }

(** User-supplied functions for the BBD preconditioner.

    @idas <node5#sss:idabbdpre> IDABBDLocalFn
    @idas <node5#sss:idabbdpre> IDABBDCommFn *)
type callbacks =
  {
    local_fn : float -> data -> data -> data -> unit;
      (** [gloc t y gd] computes [g(t, y)] into [gd]. This function
          should raise {!Sundials.RecoverableFailure} on a recoverable error,
          any other exception is treated as an unrecoverable error. *)

    comm_fn  : (float -> data -> data -> unit) option;
      (** [cfn t y] performs all interprocess communication necessary
          for the execution of [local_fn] using the input vector [y]. This
          function should raise {!Sundials.RecoverableFailure} on a recoverable
          error, any other exception is treated as an unrecoverable error. *)
  }

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
val prec_left : ?dqrely:float -> bandwidths -> callbacks
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

