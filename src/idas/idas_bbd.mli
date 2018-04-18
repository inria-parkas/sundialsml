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

(** Parallel band-block-diagonal preconditioners for IDAS (requires MPI).

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)
    @idas <node7#SECTION00742000000000000000> Using the band-block-diagonal preconditioner IDABBDPRE
 *)

(** Alias for sessions based on parallel nvectors. *)
type parallel_bsession =
  (Nvector_parallel.data, Nvector_parallel.kind) Idas.Adjoint.bsession

(** Alias for preconditioners based on parallel nvectors. *)
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

(** Functions that approximate backward residual functions using local
    computations only.  They are passed the arguments:
    - [args], the current values of forward and backward variables
              and their derivatives, and
    - [resb], a vector for storing a local approximation to the backward
              residual function {% $F_B(t, y, \dot{y}, y_B, \dot{y}_B)$ %}.

    The approximation is allowed to coincide with the actual backward
    residual function.

    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @idas <node7#SECTION00742200000000000000> IDABBDLocalFnB *)
type local_fn = Nvector_parallel.data Idas.Adjoint.bresfn_args
                -> Nvector_parallel.data
                -> unit

(** Functions that perform the interprocess communication necessary
    for the execution of {!local_fn}.
    - [args], the current values of forward and backward variables, and,
    - [resb], a vector for storing a local approximation to the backward
              function {% $F_B(t, y, \dot{y}, y_B, \dot{y}_B)$ %}.

    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @idas <node7#SECTION00742200000000000000> IDABBDCommFnB *)
type comm_fn = Nvector_parallel.data Idas.Adjoint.bresfn_args
               -> unit

(** Left preconditioning using the Parallel Band-Block-Diagonal
    module.  The difference quotient operation is controlled by
    [?dqrely], which specifies the relative increment in components of
    [y], and {!bandwidths}.

    @idas <node7#sss:lin_solv_b> IDASpgmrB
    @idas <node7#SECTION00742100000000000000> IDABBDPrecInitB *)
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

    @idas <node7#SECTION00742000000000000000> IDABBDPrecReInitB *)
val reinit : parallel_bsession -> ?dqrely:float -> int -> int -> unit

(** Returns the sizes of the real and integer workspaces used by the
    BBD preconditioner.

    @ida <node7#SECTION00742000000000000000> IDABBDPrecGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space : parallel_bsession -> int * int

(** Returns the number of calls to the right-hand side function due to
    finite difference banded Jacobian approximation in the setup function.

    @ida <node7#SECTION00742000000000000000> IDABBDPrecGetNumGfnEvals *)
val get_num_gfn_evals : parallel_bsession -> int

