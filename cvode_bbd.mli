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

(** Interface to the CVODE Parallel Band-Block-Diagonal preconditioners.

    @version VERSION()
    @author Timothy Bourke (Inria)
    @author Jun Inoue (Inria)
    @author Marc Pouzet (LIENS)
    @cvode <node5#sss:cvbbdpre> Parallel band-block-diagonal preconditioner module *)

type data = Nvector_parallel.data
type parallel_session = (data, Nvector_parallel.kind) Cvode.session
type parallel_linear_solver = (data, Nvector_parallel.kind) Cvode.linear_solver

(** TODO:
    4. Implement examples (cvode/parallel and cvodes/parallel).
    5. Update the TODO file.
  *)


type bandwidths =
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

    @cvodes <node5#sss:cvbbdpre> CVBBDLocalFn
    @cvodes <node5#sss:cvbbdpre> CVBBDCommFn *)
type callbacks =
  {
    local_fn : float -> data -> data -> unit;
      (** [gloc t y gd] computes [g(t, y)] into [gd]. This function
          should raise {!Sundials.RecoverableFailure} on a recoverable error,
          any other exception is treated as an unrecoverable error. *)

    comm_fn  : (float -> data -> unit) option;
      (** [cfn t y] performs all interprocess communication necessary
          for the execution of [local_fn] using the input vector [y]. This
          function should raise {!Sundials.RecoverableFailure} on a recoverable
          error, any other exception is treated as an unrecoverable error. *)
  }

(** Same as {!Cvode.Spils.spgmr} but with the Parallel Band-Block-Diagonal
    preconditioner.
    
    The arguments specify the maximum dimension of the Krylov subspace (pass
    [None] to use the default value [5].), the preconditioning type, the
    bandwidths described under {!bandwidths}, the relative increment in
    components of [y] used in the difference quotient approximations (pass
    [None] to use the default value [sqrt unit_roundoff]), and the callbacks
    described under {!callbacks}.

    @cvode <node5#sss:lin_solve_init> CVSpgmr
    @cvode <node5#sss:cvbbdpre> CVBBDPrecInit *)
val spgmr : int option -> Spils.preconditioning_type -> bandwidths
                -> float option -> callbacks -> parallel_linear_solver

(** Same as {!Cvode.spils.spbcg} but with the Parallel Band-Block-Diagonal
    preconditioner. The arguments are the same as for [spgmr].

    @cvode <node5#sss:lin_solve_init> CVSpbcg
    @cvode <node5#sss:cvbbdpre> CVBBDPrecInit *)
val spbcg : int option -> Spils.preconditioning_type -> bandwidths
                -> float option -> callbacks -> parallel_linear_solver

(** Same as {!Cvode.spils.spbcg} but with the Parallel Band-Block-Diagonal
    preconditioner. The arguments are the same as for [spgmr].

    @cvode <node5#sss:lin_solve_init> CVSptfqmr
    @cvode <node5#sss:cvbbdpre> CVBBDPrecInit *)
val sptfqmr : int option -> Spils.preconditioning_type -> bandwidths
                -> float option -> callbacks -> parallel_linear_solver

(** [reinit s mudq mldq dqrely] reinitializes the BBD preconditioner
    with upper ([mudq]) and lower ([mldq]) half-bandwidths to be used in the
    difference quotient Jacobian approximation, and an optional relative
    increment in components of [y] (passing [None] uses the default value [sqrt
    unit_roundoff]).

    @cvode <node5#sss:cvbbdpre> CVBBDPrecReInit *)
val reinit : parallel_session -> int -> int -> float option -> unit

(** {4 Optional output functions} *)

(** Returns the sizes of the real and integer workspaces used by the
    band-block-diagonal preconditioner module.

    @cvode <node5#sss:cvbbdpre> CVBBDPrecGetWorkSpace
    @return ([real_size], [integer_size]) *)
val get_work_space : parallel_session -> int * int

(** Returns the number of calls made to the user-supplied right-hand
    side function due to finite difference banded Jacobian approximation in the
    preconditioner setup function.

    @cvode <node5#sss:cvbbdpre> CVBBDPrecGetNumGfnEvals *)
val get_num_gfn_evals : parallel_session -> int

