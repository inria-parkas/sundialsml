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

(** Interface to the IDA Parallel Band-Block-Diagonal preconditioners.

    @version VERSION()
    @author Timothy Bourke (Inria)
    @author Jun Inoue (Inria)
    @author Marc Pouzet (LIENS)
    @ida <node5#sss:idabbdpre> Parallel band-block-diagonal preconditioner module *)

type data = Nvector_parallel.data
type kind = Nvector_parallel.kind
type parallel_session = (data, kind) Ida.session
type parallel_linear_solver = (data, kind) Ida.linear_solver

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

(** Same as {!Ida.Spils.spgmr} but with the Parallel
    Band-Block-Diagonal preconditioner.  Called like [spgmr ~maxl:maxl
    ~max_restarts:maxr ~dqrely:dqrely bandwidths cbs], where:

    - [~maxl] is the maximum dimension of the Krylov subspace.
      Defaults to [5].
    - [~max_restarts] is the maximum number of restarts.  Defaults to [5].
      Passing [0] disables restarts.
    - [~dqrely] is the relative increment in components of [y] used in
      the difference quotient approximations.  Defaults to [sqrt
      unit_roundoff].
    - [bandwidths] give the bandwidths of the preconditioning matrix.
    - [cbs] is a set of {!callbacks}.

    @ida <node5#sss:lin_solv_init> IDASpgmr
    @ida <node5#sss:idabbdpre> IDABBDPrecInit *)
val spgmr : ?maxl:int -> ?max_restarts:int -> ?dqrely:float
            -> bandwidths -> callbacks -> parallel_linear_solver

(** Same as {!Ida.spils.spbcg} but with the Parallel
    Band-Block-Diagonal preconditioner.  The arguments are the same as
    for {!spgmr}, except the maximum number of restarts
    ([~max_restarts]) cannot be specified.

    @ida <node5#sss:lin_solv_init> IDASpbcg
    @ida <node5#sss:idabbdpre> IDABBDPrecInit *)
val spbcg : ?maxl:int -> ?dqrely:float
            -> bandwidths -> callbacks -> parallel_linear_solver

(** Same as {!Ida.spils.spbcg} but with the Parallel
    Band-Block-Diagonal preconditioner.  The arguments are the same as
    for {!spgmr}, except the maximum number of restarts
    ([~max_restarts]) cannot be specified.

    @ida <node5#sss:lin_solv_init> IDASptfqmr
    @ida <node5#sss:idabbdpre> IDABBDPrecInit *)
val sptfqmr : ?maxl:int -> ?dqrely:float
            -> bandwidths -> callbacks -> parallel_linear_solver

(** [reinit s mudq mldq dqrely] reinitializes the BBD preconditioner
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

