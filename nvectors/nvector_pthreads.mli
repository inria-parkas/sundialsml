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

(** The Pthreads nvectors of Sundials (requires pthreads). 

    @version VERSION()
    @author Timothy Bourke (Inria/ENS)
    @author Jun Inoue (Inria/ENS)
    @author Marc Pouzet (UPMC/ENS/Inria)
    @cvode <node7#ss:nvec_par> NVECTOR_PARALLEL *)

(** Pthreads nvectors are based on
    {{:OCAML_DOC_ROOT(Bigarray.Array1.html)} bigarrays} of floats. *)
type data = Sundials.RealArray.t

(** Pthreads nvectors can be used wherever Serial nvectors can. *)
type kind = Nvector_serial.kind

(** The type of Pthreads nvectors. *)
type t = (data, kind) Nvector.t

(** [make nthreads n iv] creates a new Pthreads nvector with [nthreads]
    threads and [n] elements inialized to [iv]. *)
val make : int -> int -> float -> t

(** [wrap nthreads a] creates a new Pthreads nvector with [nthreads] threads
    over the elements of [a]. *)
val wrap : int -> Sundials.RealArray.t -> t

(** Aliases {!Nvector.unwrap}. *)
val unwrap : t -> Sundials.RealArray.t

(** Returns the number of threads used within a Pthreads nvector. *)
val num_threads : t -> int

(** Underlyling nvector operations on Pthreads nvectors. *)
module Ops : Nvector.NVECTOR_OPS with type t = t

