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

(** Interface to the KINSOL Parallel Band-Block-Diagonal preconditioners.

    @version VERSION()
    @author Timothy Bourke (Inria)
    @author Jun Inoue (Inria)
    @author Marc Pouzet (LIENS)
    @kinsol <node5#sss:kinbbdpre> Parallel band-block-diagonal preconditioner module *)

type data = Nvector_parallel.data
type kind = Nvector_parallel.kind
type parallel_session = (data, kind) Kinsol.session
type parallel_linear_solver = (data, kind) Kinsol_impl.linear_solver

type bandwidths =
  {
    mudq    : int; (** Upper half-bandwidth to be used in the difference
                       quotient Jacobian approximation. *)
    mldq    : int; (** Lower half-bandwidth to be used in the difference
                       quotient Jacobian approximation. *)
    mukeep  : int; (** Upper half-bandwidth of the retained banded
                       approximate Jacobian block. *)
    mlkeep  : int; (** Lower half-bandwidth of the retained banded
                       approximate Jacobian block. *)
  }

(** User-supplied functions for the BBD preconditioner.

    @kinsol <node5#sss:kinbbdpre> KINLocalFn
    @kinsol <node5#sss:kinbbdpre> KINCommFn *)
type callbacks =
  {
    local_fn : data -> data -> unit;
      (** [gloc u gval] computes [g(u)] into [gval]. This function
          should raise {!Sundials.RecoverableFailure} on a recoverable error,
          any other exception is treated as an unrecoverable error. *)

    comm_fn  : (data -> unit) option;
      (** [cfn u] performs all interprocess communication necessary
          for the execution of [local_fn] using the input vector [u]. This
          function should raise {!Sundials.RecoverableFailure} on a recoverable
          error, any other exception is treated as an unrecoverable error. *)
  }

(** Same as {!Kinsol.Spils.spgmr} but with the Parallel Band-Block-Diagonal
    preconditioner.

    Krylov iterative solver with the scaled preconditioned GMRES method. The
    arguments specify the maximum dimension of the Krylov subspace (Pass None to
    use the default value 5), the maximum number of restarts of the iterative
    linear solver (Pass None to use the default value 0), the bandwidths
    described above, the relative increment in components of [y]
    used in the difference quotient approximations (pass [None] to use the
    default value [sqrt unit_roundoff]), and the callbacks described under
    {!callbacks}.

    @kinsol <node5#sss:lin_solv_init> KINSpgmr
    @kinsol <node5#sss:kinbbdpre> KINBBDPrecInit *)
val spgmr : ?maxl:int -> ?max_restarts:int -> ?dqrely:float
          -> bandwidths -> callbacks -> parallel_linear_solver

(** Same as {!Kinsol.Spils.spbcg} but with the Parallel Band-Block-Diagonal
    preconditioner. The arguments are the same as for {!spgmr} except that a
    maximum number of restarts is not given.

    @kinsol <node5#sss:lin_solv_init> KINSpbcg
    @kinsol <node5#sss:kinbbdpre> KINBBDPrecInit *)
val spbcg : ?maxl:int -> ?dqrely:float
          -> bandwidths -> callbacks -> parallel_linear_solver

(** Same as {!Kinsol.Spils.sptfqmr} but with the Parallel Band-Block-Diagonal
    preconditioner. The arguments are the same as for {!spgmr} except that a
    maximum number of restarts is not given.

    @kinsol <node5#sss:lin_solv_init> KINSptfqmr
    @kinsol <node5#sss:kinbbdpre> KINBBDPrecInit *)
val sptfqmr : ?maxl:int -> ?dqrely:float
            -> bandwidths -> callbacks -> parallel_linear_solver

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

