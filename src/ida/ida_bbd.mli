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

(** Parallel band-block-diagonal preconditioners for IDA (requires MPI).

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)
    @ida <node5#sss:idabbdpre> Parallel band-block-diagonal preconditioner module *)

open Sundials

(** Alias for sessions based on parallel nvectors. *)
type parallel_session =
  (Nvector_parallel.data, Nvector_parallel.kind) Ida.session

(** Alias for preconditioners based on parallel nvectors. *)
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

(** Approximates the right-hand side function using local computations only.
    In the call [gloc t y y' r], [t] is the independent variable (time),
    [y] is the dependent variable vector, [y'] is the derivative vector,
    and [r] is the output vector to fill with the computed residual.

    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @idas <node5#sss:idabbdpre> IDABBDLocalFn *)
type local_fn = float
                -> Nvector_parallel.data
                -> Nvector_parallel.data
                -> Nvector_parallel.data
                -> unit

(** Functions that perform the interprocess communication necessary
    for the execution of {!local_fn}.
    In the call [cfn t y y'], [t] is the independent variable (time), [y] is
    the dependent variable vector, and [y'] is the derivative vector.

    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @idas <node5#sss:idabbdpre> IDABBDCommFn *)
type comm_fn = float -> Nvector_parallel.data -> Nvector_parallel.data -> unit

(** Left preconditioning using the Parallel Band-Block-Diagonal
    module.  The difference quotient operation is controlled by
    [?dqrely], which specifies the relative increment in components of
    [y], and {!bandwidths}.

    @ida <node5#sss:lin_solv_init> IDASpgmr
    @ida <node5#sss:idabbdpre> IDABBDPrecInit *)
val prec_left : ?dqrely:float
                -> bandwidths
                -> ?comm:comm_fn
                -> local_fn
                -> parallel_preconditioner

(** Reinitializes some BBD preconditioner parameters.
    In the call, [reinit s ~dqrely:dqrely mudq mldq], [dqrely] is the relative
    increment in the components of [y], and [mudq] and [mldq] are, respectively,
    the upper-half and lower-half bandwidths of the difference quotient
    Jacobian approximation.

    @ida <node5#sss:idabbdpre> IDABBDPrecReInit *)
val reinit : parallel_session -> ?dqrely:float -> int -> int -> unit

(** Returns the sizes of the real and integer workspaces used by the
    BBD preconditioner.

    @ida <node5#sss:idabbdpre> IDABBDPrecGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space : parallel_session -> int * int

(** Returns the number of calls to the right-hand side function due to
    finite difference banded Jacobian approximation in the setup function.

    @ida <node5#sss:idabbdpre> IDABBDPrecGetNumGfnEvals *)
val get_num_gfn_evals : parallel_session -> int

