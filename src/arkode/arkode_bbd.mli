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
    @arkode <Usage/ARKStep_c_interface/Preconditioners.html#a-parallel-band-block-diagonal-preconditioner-module> Parallel band-block-diagonal preconditioner module *)

(** Alias for sessions based on parallel nvectors. *)
type 'step parallel_session =
  (Nvector_parallel.data, Nvector_parallel.kind, 'step) Arkode_impl.session
  constraint 'step = [<Arkode.arkstep|Arkode.mristep]

(** Alias for preconditioners based on parallel nvectors. *)
type 'step parallel_preconditioner =
  (Nvector_parallel.data, Nvector_parallel.kind, 'step) Arkode.Spils.preconditioner
  constraint 'step = [<Arkode.arkstep|Arkode.mristep]

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

    @arkode_precond ARKLocalFn *)
type local_fn = float -> Nvector_parallel.data -> Nvector_parallel.data -> unit

(** Functions that perform the interprocess communication necessary
    for the execution of {!local_fn}.
    In the call [cfn t y], [t] is the independent variable (time) and [y] is
    the input vector.

    Raising {!Sundials.RecoverableFailure} signals a recoverable error.
    Other exceptions signal unrecoverable errors.

    @arkode_precond ARKCommFn *)
type comm_fn = float -> Nvector_parallel.data -> unit

(** Left preconditioning using the Parallel Band-Block-Diagonal
    module.  The difference quotient operation is controlled by
    [?dqrely], which specifies the relative increment in components of
    [y], and {!bandwidths}.

    NB: Band-Block-Diagonal preconditioners may not be used for problems
    involving a non-identity mass matrix.

    @arkode_precond ARKBBDPrecInit *)
val prec_left : ?dqrely:float
                -> bandwidths
                -> ?comm:comm_fn
                -> local_fn
                -> 's parallel_preconditioner

(** Right preconditioning using the Parallel Band-Block-Diagonal
    module.  The difference quotient operation is controlled by
    [?dqrely], which specifies the relative increment in components of
    [y], and {!bandwidths}.

    @arkode_precond ARKBBDPrecInit *)
val prec_right : ?dqrely:float
                 -> bandwidths
                 -> ?comm:comm_fn
                 -> local_fn
                 -> 's parallel_preconditioner

(** Preconditioning from both sides using the Parallel Band-Block-Diagonal
    module.
    The difference quotient operation is controlled by [?dqrely],
    the relative increment in components of [y], and {!bandwidths}.

    @arkode_precond ARKBBDPrecInit *)
val prec_both : ?dqrely:float
                -> bandwidths
                 -> ?comm:comm_fn
                 -> local_fn
                -> 's parallel_preconditioner

(** Reinitializes some BBD preconditioner parameters.
    In the call, [reinit s ~dqrely:dqrely mudq mldq], [dqrely] is the relative
    increment in the components of [y], and [mudq] and [mldq] are, respectively,
    the upper-half and lower-half bandwidths of the difference quotient
    Jacobian approximation.

    @arkode_precond ARKBBDPrecReInit *)
val reinit : 's parallel_session -> ?dqrely:float -> int -> int -> unit

(** Returns the sizes of the real and integer workspaces used by the
    BBD preconditioner.

    @arkode_precond ARKBBDPrecGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space : 's parallel_session -> int * int

(** Returns the number of calls to the right-hand side function due to
    finite difference banded Jacobian approximation in the setup function.

    @arkode_precond ARKBBDPrecGetNumGfnEvals *)
val get_num_gfn_evals : 's parallel_session -> int

