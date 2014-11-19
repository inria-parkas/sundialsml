(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(** Parallel band-block-diagonal preconditioners for KINSOL (requires MPI).

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)
    @kinsol <node5#sss:kinbbdpre> Parallel band-block-diagonal preconditioner module *)

(** Alias for sessions based on parallel nvectors. *)
type parallel_session =
      (Nvector_parallel.data, Nvector_parallel.kind) Kinsol.session

(** Alias for preconditioners based on parallel nvectors. *)
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

(** Approximates the system function using local computations only.
    In the call [gloc u gval], [u] is the current estimated solution
    and [gval] stores the computed solution.

    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @kinsol <node5#sss:kinbbdpre> KINLocalFn *)
type local_fn = Nvector_parallel.data -> Nvector_parallel.data -> unit

(** Functions that perform the interprocess communication necessary
    for the execution of {!local_fn}. In the call [cfn u], [u] is the input
    vector.

    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @kinsol <node5#sss:kinbbdpre> KINCommFn *)
type comm_fn = Nvector_parallel.data -> unit

(** Right preconditioning using the Parallel Band-Block-Diagonal module.
    The difference quotient operation is controlled by [?dqrely],
    the relative increment in components of [y], and {!bandwidths}.

    @kinsol <node5#sss:lin_solv_init> KINSpgmr
    @kinsol <node5#sss:kinbbdpre> KINBBDPrecInit *)
val prec_right : ?dqrely:float
                 -> bandwidths
                 -> ?comm_fn:comm_fn
                 -> local_fn
                 -> parallel_preconditioner

(** Returns the sizes of the real and integer workspaces used by the
    BBD preconditioner.

    @kinsol <node5#sss:kinbbdpre> KINBBDPrecGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space : parallel_session -> int * int

(** Returns the number of calls to the right-hand side function due to
    finite difference banded Jacobian approximation in the setup function.

    @kinsol <node5#sss:kinbbdpre> KINBBDPrecGetNumGfnEvals *)
val get_num_gfn_evals : parallel_session -> int

