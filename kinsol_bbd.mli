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

(** Parallel band-block-diagonal preconditioners for KINSOL (requires MPI).

    @version VERSION()
    @author Timothy Bourke (Inria)
    @author Jun Inoue (Inria)
    @author Marc Pouzet (LIENS)
    @kinsol <node5#sss:kinbbdpre> Parallel band-block-diagonal preconditioner module *)

(** An alias for sessions based on parallel nvectors. *)
type parallel_session =
      (Nvector_parallel.data, Nvector_parallel.kind) Kinsol.session

(** An alias for preconditioners based on parallel nvectors. *)
type parallel_preconditioner =
      (Nvector_parallel.data, Nvector_parallel.kind) Kinsol.Spils.preconditioner

(** The bandwidths for the difference quotient Jacobian operation. *)
type bandwidths =
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

(** [gloc u gval] computes [g(u)] into [gval].
    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @kinsol <node5#sss:kinbbdpre> KINLocalFn *)
type local_fn = Nvector_parallel.data -> Nvector_parallel.data -> unit

(** [cfn u] performs all interprocess communication necessary
    for the execution of [local_fn] using the input vector [u].
    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @kinsol <node5#sss:kinbbdpre> KINCommFn *)
type comm_fn = Nvector_parallel.data -> unit

(** Same as {!Kinsol.Spils.prec_right} but sets up the Parallel
    Band-Block-Diagonal preconditioner included in KINSOL.  Called
    like [prec_right ~dqrely:dqrely callbacks], where:

    - [~dqrely] gives the relative increment in components of [y] used in
      the difference quotient approximations
      (defaults to [sqrt unit_roundoff]).
    - [bandwidths] specify the bandwidths to be used in the difference
      quotient Jacobian operation.
    - [callbacks] gives the preconditioning callbacks.  See the
      {!callbacks} type.

    @kinsol <node5#sss:lin_solv_init> KINSpgmr
    @kinsol <node5#sss:kinbbdpre> KINBBDPrecInit *)
val prec_right : ?dqrely:float
                 -> bandwidths
                 -> ?comm_fn:comm_fn
                 -> local_fn
                 -> parallel_preconditioner

(** {4 Optional output functions} *)

(** Returns the sizes of the real and integer workspaces used by the
    band-block-diagonal preconditioner module.

    @kinsol <node5#sss:kinbbdpre> KINBBDPrecGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space : parallel_session -> int * int

(** Returns the number of calls made to the user-supplied right-hand
    side function due to finite difference banded Jacobian approximation in
    the preconditioner setup function.

    @kinsol <node5#sss:kinbbdpre> KINBBDPrecGetNumGfnEvals *)
val get_num_gfn_evals : parallel_session -> int

