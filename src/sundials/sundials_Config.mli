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

(** Information on Sundials/ML and the underlying Sundials library. *)

(** The [major], [minor], [patch], and [binding] version numbers of
    Sundials/ML.
    The first three elements correspond to the maximum supported version
    of the underlying Sundials/C library.
    The [binding] number distinguishes updates to the binding (restarting
    from 0 for each increment of the other elements). *)
val version : int * int * int * int

(** The [major], [minor], and [patch] version numbers of the underlying
    Sundials/C library. The OCaml interface may have been built against
    an older version of the Sundials/C library  *)
val sundials_version : int * int * int

(** Indicates whether the interface was compiled with BLAS/LAPACK support. *)
val lapack_enabled : bool

(** Indicates whether the underlying library was built with support for
    monitoring functions. *)
val monitoring_enabled : bool

(** Indicates whether the parallel nvectors and linear solvers are available. *)
val mpi_enabled : bool

(** Indicates whether the KLU sparse linear solver is available. *)
val klu_enabled : bool

(** Indicates whether the SuperLU_MT sparse linear solver is available. *)
val superlumt_enabled : bool

(** Indicates whether pthreads-based nvectors are available. *)
val nvecpthreads_enabled : bool

(** Indicates whether openmp-based nvectors are available. *)
val nvecopenmp_enabled : bool

(** The largest value representable as a real.

    @cvode <node5#s:types> Data Types *)
val big_real : float

(** The smallest value representable as a real.

    @cvode <node5#s:types> Data Types *)
val small_real : float

(** The difference bewteen 1.0 and the minimum real greater than 1.0.

    @cvode <node5#s:types> Data Types *)
val unit_roundoff : float

(** Raised for features that are not available in the installed version of
    the underlying sundials library. See {!sundials_version}. *)
exception NotImplementedBySundialsVersion

