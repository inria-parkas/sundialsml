(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2015 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(** Parallel band-block-diagonal preconditioners for ARKODE (requires MPI).

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)
    @noarkode <node5#sss:arkbbdpre> Parallel band-block-diagonal preconditioner module *)

(** Alias for sessions based on parallel nvectors. *)
type parallel_session =
      (Nvector_parallel.data, Nvector_parallel.kind) Arkode.session

(** Alias for preconditioners based on parallel nvectors. *)
type parallel_preconditioner =
      (Nvector_parallel.data, Nvector_parallel.kind) Arkode.Spils.preconditioner

(** The bandwidths for the difference quotient Jacobian operation. *)
type bandwidths = Arkode_impl.ArkodeBbdTypes.bandwidths =
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
    That is, it computes the approximation {% $g(t,y) \approx f_I(t,y)$%}.
    In the call [gloc t y g], [t] is the independent variable (time), [y] is
    the input vector, and [g] stores the computed derivatives.

    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @noarkode <node5#sss:arkbbdpre> ARKLocalFn *)
type local_fn = float -> Nvector_parallel.data -> Nvector_parallel.data -> unit

(** Functions that perform the interprocess communication necessary
    for the execution of {!local_fn}.
    In the call [cfn t y], [t] is the independent variable (time) and [y] is
    the input vector.

    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @noarkode <node5#sss:arkbbdpre> ARKCommFn *)
type comm_fn = float -> Nvector_parallel.data -> unit

(** Left preconditioning using the Parallel Band-Block-Diagonal
    module.  The difference quotient operation is controlled by
    [?dqrely], which specifies the relative increment in components of
    [y], and {!bandwidths}.

    @noarkode <node5#sss:arkbbdpre> ARKBBDPrecInit *)
val prec_left : ?dqrely:float
                -> bandwidths
                -> ?comm_fn:comm_fn
                -> local_fn
                -> parallel_preconditioner

(** Right preconditioning using the Parallel Band-Block-Diagonal
    module.  The difference quotient operation is controlled by
    [?dqrely], which specifies the relative increment in components of
    [y], and {!bandwidths}.

    @noarkode <node5#sss:arkbbdpre> ARKBBDPrecInit *)
val prec_right : ?dqrely:float
                 -> bandwidths
                 -> ?comm_fn:comm_fn
                 -> local_fn
                 -> parallel_preconditioner

(** Preconditioning from both sides using the Parallel Band-Block-Diagonal
    module.
    The difference quotient operation is controlled by [?dqrely],
    the relative increment in components of [y], and {!bandwidths}.

    @noarkode <node5#sss:arkbbdpre> ARKBBDPrecInit *)
val prec_both : ?dqrely:float
                -> bandwidths
                 -> ?comm_fn:comm_fn
                 -> local_fn
                -> parallel_preconditioner

(** Reinitializes some BBD preconditioner parameters.
    In the call, [reinit s ~dqrely:dqrely mudq mldq], [dqrely] is the relative
    increment in the components of [y], and [mudq] and [mldq] are, respectively,
    the upper-half and lower-half bandwidths of the difference quotient
    Jacobian approximation.

    @noarkode <node5#sss:arkbbdpre> ARKBBDPrecReInit *)
val reinit : parallel_session -> ?dqrely:float -> int -> int -> unit

(** Returns the sizes of the real and integer workspaces used by the
    BBD preconditioner.

    @noarkode <node5#sss:arkbbdpre> ARKBBDPrecGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space : parallel_session -> int * int

(** Returns the number of calls to the right-hand side function due to
    finite difference banded Jacobian approximation in the setup function.

    @noarkode <node5#sss:arkbbdpre> ARKBBDPrecGetNumGfnEvals *)
val get_num_gfn_evals : parallel_session -> int

