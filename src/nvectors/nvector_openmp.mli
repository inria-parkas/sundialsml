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

(** The OpenMP nvectors of Sundials (requires OpenMP). 

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)
    @cvode <node7#ss:nvec_par> NVECTOR_PARALLEL *)

(** OpenMP nvectors are based on
    {{:OCAML_DOC_ROOT(Bigarray.Array1.html)} bigarrays} of floats. *)
type data = Sundials.RealArray.t

(** Represents the internal layout of an OpenMP nvector.
    OpenMP nvectors can usually be used wherever serial nvectors can. *)
type kind = [`OpenMP|Nvector_serial.kind]

(** The type of OpenMP nvectors. *)
type t = (data, kind) Nvector.t

(** [make nthreads n iv] creates a new OpenMP nvector with [nthreads]
    threads and [n] elements inialized to [iv]. *)
val make : int -> int -> float -> t

(** [wrap nthreads a] creates a new OpenMP nvector with [nthreads] threads
    over the elements of [a]. *)
val wrap : int -> Sundials.RealArray.t -> t

(** Aliases {!Nvector.unwrap}. *)
val unwrap : t -> Sundials.RealArray.t

(** Pretty-print an OpenMP nvector using the
    {{:OCAML_DOC_ROOT(Format.html)} Format} module. *)
val pp : Format.formatter -> t -> unit

(** Returns the number of threads used within an OpenMP nvector. *)
val num_threads : t -> int

(** Underlying nvector operations on OpenMP nvectors. *)
module Ops : Nvector.NVECTOR_OPS with type t = t

