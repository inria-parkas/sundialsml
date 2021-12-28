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

(** Parallel band-block-diagonal preconditioners for CVODES (requires MPI).

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)
    @cvodes <Usage/ADJ.html?highlight=CVBBDPrecInitB#using-the-band-block-diagonal-preconditioner-cvbbdpre> Using the band-block-diagonal preconditioner CVBBDPRE
 *)

(** Alias for sessions based on parallel nvectors. *)
type parallel_bsession =
      (Nvector_parallel.data, Nvector_parallel.kind) Cvodes.Adjoint.bsession

(** Alias for preconditioners based on parallel nvectors. *)
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

(** Functions that approximate backward right-hand side functions
    using local computations only.  They are passed the arguments:
    - [args], the current values of forward and backward variables, and,
    - [yb'], a vector for storing a local approximation to the backward
             right-hand side function {% $\dot{y}_B = f_B(t, y, y_B)$%}.

    The approximation is allowed to coincide with the actual backward
    right-hand side function.

    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @cvodes_adj CVBBDLocalFnB *)
type local_fn = Nvector_parallel.data Cvodes.Adjoint.brhsfn_args
                -> Nvector_parallel.data
                -> unit

(** Functions that perform the interprocess communication necessary
    for the execution of {!local_fn}.  They are passed the current
    values of forward and backward variables.

    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @cvodes_adj CVBBDCommFnB *)
type comm_fn = Nvector_parallel.data Cvodes.Adjoint.brhsfn_args -> unit

(** Left preconditioning using the Parallel Band-Block-Diagonal
    module.  The difference quotient operation is controlled by
    [?dqrely], which specifies the relative increment in components of
    [y], and {!bandwidths}.

    @cvodes_adj CVBBDPrecInitB *)
val prec_left : ?dqrely:float
                -> bandwidths
                -> ?comm:comm_fn
                -> local_fn
                -> parallel_preconditioner

(** Right preconditioning using the Parallel Band-Block-Diagonal
    module.  The difference quotient operation is controlled by
    [?dqrely], which specifies the relative increment in components of
    [y], and {!bandwidths}.

    @cvodes_adj CVBBDPrecInitB *)
val prec_right : ?dqrely:float
                 -> bandwidths
                 -> ?comm:comm_fn
                 -> local_fn
                 -> parallel_preconditioner

(** Preconditioning from both sides using the Parallel Band-Block-Diagonal
    module.
    The difference quotient operation is controlled by [?dqrely],
    the relative increment in components of [y], and {!bandwidths}.

    @cvodes_adj CVBBDPrecInitB *)
val prec_both : ?dqrely:float
                -> bandwidths
                -> ?comm:comm_fn
                -> local_fn
                -> parallel_preconditioner

(** Reinitializes some BBD preconditioner parameters.
    In the call, [reinit s ~dqrely:dqrely mudq mldq], [dqrely] is the relative
    increment in the components of [y], and [mudq] and [mldq] are, respectively,
    the upper-half and lower-half bandwidths of the difference quotient
    Jacobian approximation.

    @cvodes_adj CVBBDPrecReInitB *)
val reinit : parallel_bsession -> ?dqrely:float -> int -> int -> unit

(** Returns the sizes of the real and integer workspaces used by the
    BBD preconditioner.

    @cvodes CVBBDPrecGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space : parallel_bsession -> int * int

(** Returns the number of calls to the right-hand side function due to
    finite difference banded Jacobian approximation in the setup function.

    @cvodes CVBBDPrecGetNumGfnEvals *)
val get_num_gfn_evals : parallel_bsession -> int

