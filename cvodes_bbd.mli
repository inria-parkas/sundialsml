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

(** Interface to the CVODES Adjoint Parallel Band-Block-Diagonal preconditioner.

    @version VERSION()
    @author Timothy Bourke (Inria)
    @author Jun Inoue (Inria)
    @author Marc Pouzet (LIENS)
    @cvodes <node7#SECTION00742000000000000000> Using the band-block-diagonal preconditioner CVBBDPRE
 *)
type data = Nvector_parallel.data
type kind = Nvector_parallel.kind
type parallel_bsession = (data, kind) Cvodes.Adjoint.bsession
type parallel_linear_solver = (data, kind) Cvodes.Adjoint.linear_solver

type bandwidths = Cvode_bbd.bandwidths =
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

    @cvodes <node7#SECTION00742200000000000000> CVBBDLocalFnB
    @cvodes <node7#SECTION00742200000000000000> CVBBDCommFnB *)
type callbacks =
  {
    local_fn : float -> data -> data -> data -> unit;
      (** [gloc t y yb gb] computes [g(t, y)] into [gb] from the value of the
          independent variable [t], the current forward solution vector [y], and
          the current value of the backward dependent variable vector. This
          function should raise {!Sundials.RecoverableFailure} on a recoverable
          error, any other exception is treated as an unrecoverable error. *)

    comm_fn  : (float -> data -> data -> unit) option;

      (** [cfn t y yb] performs all interprocess communication necessary for the
          execution of [local_fn] using the forward solution vector [y] and the
          backward dependent variable vector [yb]. This function should raise
          {!Sundials.RecoverableFailure} on a recoverable error, any other
          exception is treated as an unrecoverable error. *)
  }

(** Same as {!Cvode.Spils.spgmr} but with the Parallel
    Band-Block-Diagonal preconditioner.  Called like [spgmr ~maxl:maxl
    prec_type bandwidths ~dqrely:dqrely callbacks], where:

    - [~maxl] gives the maximum dimension of the Krylov subspace
      (defaults to [5]).
    - [prec_type] is the type of preconditioning.  It is legal to specify
      [PrecTypeNone], though in that case {!Cvodes.Spils.spgmr} would work
      just as well.
    - [bandwidths] specify the bandwidths to be used in the difference
      quotient Jacobian operation.
    - [~dqrely] gives the relative increment in components of [y] used
      in the difference quotient approximations (defaults to [sqrt
      Sundials.unit_roundoff]).
    - [callbacks] gives the preconditioning callbacks.  See the
      {!callbacks} type.

    @cvodes <node7#sss:lin_solv_b> CVSpgmrB
    @cvodes <node7#SECTION00742100000000000000> CVBBDPrecInitB *)
val spgmr : ?maxl:int -> Spils.preconditioning_type -> bandwidths
            -> ?dqrely:float -> callbacks -> parallel_linear_solver

(** Same as {!Cvodes.Adjoint.Spils.spbcg} but with the Parallel
    Band-Block-Diagonal preconditioner. The arguments are the same as for
    {!spgmr}.CVBBDPrecInit

    @cvodes <node7#sss:lin_solv_b> CVSpbcgB
    @cvodes <node7#SECTION00742100000000000000> CVBBDPrecInitB *)
val spbcg : ?maxl:int -> Spils.preconditioning_type -> bandwidths
            -> ?dqrely:float -> callbacks -> parallel_linear_solver

(** Same as {!Cvodes.Adjoint.Spils.sptfqmr} but with the Parallel
    Band-Block-Diagonal preconditioner. The arguments are the same as for
    {!spgmr}.

    @cvodes <node7#sss:lin_solv_b> CVSptfqmrB
    @cvodes <node7#SECTION00742100000000000000> CVBBDPrecInitB *)
val sptfqmr : ?maxl:int -> Spils.preconditioning_type -> bandwidths
              -> ?dqrely:float -> callbacks -> parallel_linear_solver

(** [reinit s mudq mldq ~dqrely:dqrely] reinitializes the BBD
    preconditioner with upper ([mudq]) and lower ([mldq])
    half-bandwidths to be used in the difference quotient Jacobian
    approximation, and an optional relative increment [dqrely] in
    components of [y].  [dqrely] defaults to
    [sqrt Sundials.unit_roundoff].

    @cvodes <node7#SECTION00742000000000000000> CVBBDPrecReInitB *)
val reinit : parallel_bsession -> ?dqrely:float -> int -> int -> unit

(** {4 Optional output functions} *)

(** Returns the sizes of the real and integer workspaces used by the
    band-block-diagonal preconditioner module.

    @cvode <node7#SECTION00742000000000000000> CVBBDPrecGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space : parallel_bsession -> int * int

(** Returns the number of calls made to the user-supplied right-hand
    side function due to finite difference banded Jacobian approximation in the
    preconditioner setup function.

    @cvode <node7#SECTION00742000000000000000> CVBBDPrecGetNumGfnEvals *)
val get_num_gfn_evals : parallel_bsession -> int

