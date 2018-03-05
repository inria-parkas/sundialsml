(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*             Timothy Bourke, Jun Inoue, and Marc Pouzet              *)
(*             (Inria/ENS)     (Inria/ENS)    (UPMC/ENS/Inria)         *)
(*                                                                     *)
(*  Copyright 2018 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(* TODO: Change the name/location of this module? *)

(** Generic linear solvers.

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)

    @nocvode <node> Description of the SUNLinearSolver module
    @since 3.0.0 *)

(** Direct Linear Solvers. *)
module Direct : sig
  (** A direct linear solver.

      @nocvode <node> Description of the SUNLinearSolver module
      @nocvode <node> SUNLinearSolver *)
  type ('m, 'nd, 'nk) t

  (** Alias for linear solvers that are restricted to serial nvectors. *)
  type ('m, 'nk) serial_t
    = ('m, Nvector_serial.data, [>Nvector_serial.kind] as 'nk) t

  (** Creates a direct linear solver on dense matrices. The nvector and matrix
      argument are used to determine the linear system size and to assess
      compatibility with the linear solver implementation.

    @nocvode <node> SUNDenseLinearSolver *)
  val dense :
    'k Nvector_serial.any
    -> 'k Matrix.dense
    -> (Matrix.Dense.t, 'k) serial_t

  (** Creates a direct linear solver on dense matrices using LAPACK.
      See {!make}. Only available if {!Sundials.lapack_enabled}.

    @nocvode <node> SUNLapackDense *)
  val lapack_dense :
    'k Nvector_serial.any
    -> 'k Matrix.dense
    -> (Matrix.Dense.t, 'k) serial_t

  (** Creates a direct linear solver on banded matrices. The nvector and matrix
      argument are used to determine the linear system size and to assess
      compatibility with the linear solver implementation.

    @nocvode <node> SUNBandLinearSolver *)
  val band :
    'k Nvector_serial.any
    -> 'k Matrix.band
    -> (Matrix.Band.t, 'k) serial_t

  (** Creates a direct linear solver on banded matrices using LAPACK.
      See {!make}. Only available if {!Sundials.lapack_enabled}.

    @nocvode <node> SUNLapackBand *)
  val lapack_band :
    'k Nvector_serial.any
    -> 'k Matrix.band
    -> (Matrix.Band.t, 'k) serial_t

  (** KLU direct linear solver operating on sparse matrices (requires KLU). *)
  module Klu : sig
    (** Creates a direct linear solver on sparse matrices using KLU. The
        nvector and matrix argument are used to determine the linear system
        size and to assess compatibility with the linear solver implementation.

        @raise Sundials.NotImplementedBySundialsVersion Solver not available.
        @nocvode <node> SUNKLU *)
    val make :
      'k Nvector_serial.any
      -> ('s, 'k) Matrix.sparse
      -> ('s Matrix.Sparse.t, 'k) serial_t

    (** Reinitializes memory and flags for a new factorization (symbolic and
        numeric) at the next solver setup call. In the call [reinit ls a nnz],
        [a] is the Jacobian matrix, which is reinitialized with the given
        number of non-zeros if [nnz] if given. New symbolic and numeric
        factorizations will be completed at the next solver step.

        @nocvode <node> SUNKLUReInit *)
    val reinit : ('s Matrix.Sparse.t, 'k) serial_t
      -> ('s, 'k) Matrix.sparse -> int option -> unit

    (** The ordering algorithm used for reducing fill. *)
    type ordering =
         Amd      (** Approximate minimum degree permutation. *)
       | ColAmd   (** Column approximate minimum degree permutation. *)
       | Natural  (** Natural ordering. *)

    (** Sets the ordering algorithm used to minimize fill-in.

        @nocvode <node> SUNKLUSetOrdering *)
    val set_ordering : ('s Matrix.Sparse.t, 'k) serial_t -> ordering -> unit
  end

  (** SuperLUMT direct linear solver operating on sparse matrices (requires
      SuperLUMT). *)
  module Superlumt : sig

    (** Creates a direct linear solver on sparse matrices using SuperLUMT. The
        nvector and matrix argument are used to determine the linear system
        size and to assess compatibility with the linear solver implementation.

        @raise Sundials.NotImplementedBySundialsVersion Solver not available.
        @nocvode <node> SUNSuperLUMT *)
    val make :
      'k Nvector_serial.any
      -> ('s, 'k) Matrix.sparse
      -> ('s Matrix.Sparse.t, 'k) serial_t

    (** The ordering algorithm used for reducing fill. *)
    type ordering =
         Natural       (** Natural ordering. *)
       | MinDegreeProd (** Minimal degree ordering on $J^T J$. *)
       | MinDegreeSum  (** Minimal degree ordering on $J^T + J$. *)
       | ColAmd        (** Column approximate minimum degree permutation. *)

    (** Sets the ordering algorithm used to minimize fill-in.

        @nocvode <node> SUNSuperLUMTSetOrdering *)
    val set_ordering : ('s Matrix.Sparse.t, 'k) serial_t -> ordering -> unit
  end

  (* TODO: Add custom direct lsolver (subset of generic functions). *)
end

(** Iterative linear Solvers. *)
module Iterative : sig

  (** An iterative linear solver.

      @nocvode <node> Description of the SUNLinearSolver module
      @nocvode <node> SUNLinearSolver *)
  type ('nd, 'nk, 'f) t

  (** {3:iter_solvers Solvers} *)

  (** Krylov iterative solver using the scaled preconditioned biconjugate
      stabilized (Bi-CGStab) method. The [maxl] arguments gives the maximum
      dimension of the Krylov subspace (defaults to 5). The nvector argument
      is used as a template.

      @nocvode <node> SUNSPBCGS *)
  val spbcgs : ?maxl:int -> ('d, 'k) Nvector.t -> ('d, 'k, [`Spbcgs]) t

  (** Krylov iterative solver using the scaled preconditioned flexible
      generalized minimum residual (GMRES) method. The [maxl] arguments gives
      the maximum dimension of the Krylov subspace (defaults to 5). The
      nvector argument is used as a template.

      @nocvode <node> SUNSPFGMR *)
  val spfgmr : ?maxl:int -> ?gs_type:Spils.gramschmidt_type
               -> ('d, 'k) Nvector.t -> ('d, 'k, [`Spfgmr]) t

  (** Krylov iterative solver using the scaled preconditioned generalized
      minimum residual (GMRES) method. The [maxl] arguments gives the maximum
      dimension of the Krylov subspace (defaults to 5). The nvector argument
      is used as a template.
      
      @nocvode <node> SUNSPGMR *)
  val spgmr : ?maxl:int -> ?gs_type:Spils.gramschmidt_type 
              -> ('d, 'k) Nvector.t -> ('d, 'k, [`Spgmr]) t

  (** Krylov iterative with the scaled preconditioned transpose-free
      quasi-minimal residual (SPTFQMR) method. The [maxl] arguments gives the
      maximum dimension of the Krylov subspace (defaults to 5). The nvector
      argument is used as a template.
      
      @nocvode <node> SUNSPTFQMR *)
  val sptfqmr : ?maxl:int -> ('d, 'k) Nvector.t -> ('d, 'k, [`Sptfqmr]) t

  (** Krylov iterative solver using the preconditioned conjugate gradient
      (PCG) method. The [maxl] arguments gives the maximum dimension of the
      Krylov subspace (defaults to 5). The nvector argument is used as a
      template.

      @nocvode <node> SUNPCG *)
  val pcg : ?maxl:int -> ('d, 'k) Nvector.t -> ('d, 'k, [`Pcg]) t

  (* TODO: Add custom spils lsolver (subset of generic functions). *)

  (** {3:iter_param Solver parameters} *)

  (** Updates the number of linear solver iterations to allow.

      @nocvode <node> SUNSPBCGSSetMaxl
      @nocvode <node> SUNSPTFQMRSetMaxl
      @nocvode <node> SUNPCGSetMaxl *)
  val set_maxl : ('d, 'k, [< `Spbcgs|`Sptfqmr|`Pcg]) t -> int -> unit

  (** Sets the Gram-Schmidt orthogonalization to use.

      @nocvode <node> SUNSPGMRSetGSType
      @nocvode <node> SUNSPFGMRSetGSType *)
  val set_gs_type :
    ('d, 'k, [< `Spfgmr|`Spgmr]) t -> Spils.gramschmidt_type -> unit

  (** Sets the number of GMRES restarts to allow.
 
      @nocvode <node> SUNSPGMRSetMaxRestarts
      @nocvode <node> SUNSPFGMRSetMaxRestarts *)
  val set_max_restarts : ('d, 'k, [< `Spfgmr|`Spgmr]) t -> int -> unit
end

(* TODO: Add custom linear solvers
   - custom direct lsolver (subset of generic functions)
   - custom spils lsolver (subset of generic functions)
   - custom lsolver (four functions to each solver as before)
 *)

(*
type ('a, 'd, 'k) atimes_fn
  = 'a -> ('d, 'k) Nvector.t -> ('d, 'k) Nvector.t -> unit

type ('ls, 'm, 'd, 'k) lsolver_ops = {
  m_initialize : 'ls -> unit;
  m_setup : 'ls -> 'm -> unit;
  m_solve
    : 'ls -> 'm -> ('d, 'k) Nvector.t -> ('d, 'k) Nvector.t -> float -> unit;
  m_set_atimes : ('ls -> ('a, 'd, 'k) atimes_fn -> unit) option

}
*)

