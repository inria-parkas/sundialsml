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
(*               User Documentation for IDA v1.1.0                     *)
(*                Alan C. Hindmarsh and Radu Serban                    *)
(*              Center for Applied Scientific Computing                *)
(*              Lawrence Livermore National Laboratory                 *)
(*                                                                     *)
(***********************************************************************)

(** Parallel band-block-diagonal preconditioners for IDAS (requires MPI).

    @version VERSION()
    @author Timothy Bourke (Inria)
    @author Jun Inoue (Inria)
    @author Marc Pouzet (LIENS)
    @idas <node7#SECTION00742000000000000000> Using the band-block-diagonal preconditioner IDABBDPRE
 *)

(** An alias for sessions based on parallel nvectors. *)
type parallel_bsession =
      (Nvector_parallel.data, Nvector_parallel.kind) Idas.Adjoint.bsession

(** An alias for preconditioners based on parallel nvectors. *)
type parallel_preconditioner =
      (Nvector_parallel.data, Nvector_parallel.kind)
                          Idas.Adjoint.Spils.preconditioner

(** The bandwidths for the difference quotient Jacobian operation. *)
type bandwidths = Ida_bbd.bandwidths =
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

    @idas <node7#SECTION00742200000000000000> IDABBDLocalFnB *)
type local_fn = float
                -> Nvector_parallel.data
                -> Nvector_parallel.data
                -> Nvector_parallel.data
                -> Nvector_parallel.data
                -> Nvector_parallel.data
                -> unit

(** [cfn t y yb] performs all interprocess communication necessary for the
    execution of [local_fn] using the forward solution vector [y] and the
    backward dependent variable vector [yb].
    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @idas <node7#SECTION00742200000000000000> IDABBDCommFnB *)
type comm_fn = float
               -> Nvector_parallel.data
               -> Nvector_parallel.data
               -> Nvector_parallel.data
               -> Nvector_parallel.data
               -> unit

(** Same as {!Idas.Adjoint.Spils.prec_left} but uses the Parallel
    Band-Block-Diagonal preconditioner included in IDAS.  Called like
    [prec_left ~dqrely:dqrely bandwidths callbacks], where:

    - [~dqrely] gives the relative increment in components of [y] used in
      the difference quotient approximations
      (defaults to [sqrt unit_roundoff]).
    - [bandwidths] specify the bandwidths to be used in the difference
      quotient Jacobian operation.
    - [callbacks] gives the preconditioning callbacks.  See the
      {!callbacks} type.

    @idas <node7#sss:lin_solv_b> IDASpgmrB
    @idas <node7#SECTION00742100000000000000> IDABBDPrecInitB *)
val prec_left : ?dqrely:float
                -> bandwidths
                -> ?comm_fn:comm_fn
                -> local_fn
                -> parallel_preconditioner

(** [reinit s mudq mldq dqrely] reinitializes the BBD preconditioner
    with upper ([mudq]) and lower ([mldq]) half-bandwidths to be used in the
    difference quotient Jacobian approximation, and an optional relative
    increment in components of [y] (passing [None] uses the default value [sqrt
    unit_roundoff]).

    @idas <node7#SECTION00742000000000000000> IDABBDPrecReInitB *)
val reinit : parallel_bsession -> ?dqrely:float -> int -> int -> unit

(** {4 Optional output functions} *)

(** Returns the sizes of the real and integer workspaces used by the
    band-block-diagonal preconditioner module.

    @ida <node7#SECTION00742000000000000000> IDABBDPrecGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space : parallel_bsession -> int * int

(** Returns the number of calls made to the user-supplied right-hand
    side function due to finite difference banded Jacobian approximation in the
    preconditioner setup function.

    @ida <node7#SECTION00742000000000000000> IDABBDPrecGetNumGfnEvals *)
val get_num_gfn_evals : parallel_bsession -> int

