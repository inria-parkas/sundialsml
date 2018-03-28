(***********************************************************************)
(*                                                                     *)
(*                   OCaml interface to Sundials                       *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2018 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a New BSD License, refer to the file LICENSE.                *)
(*                                                                     *)
(***********************************************************************)

(* Types underlying Lsolver.  *)

(* This module defines the lsolver types which are manipulated (abstractly)
   via the Lsolver module. They are declared here so that they can also be
   used (concretely) from the <solver> modules. To ensure that this type will
   be opaque outside of Sundials/ML, we simply do not install the
   lsolver_impl.cmi file. *)

module Klu = struct

  (* Must correspond with lsolver_ml.h:lsolver_klu_ordering_tag *)
  type ordering =
    Amd
  | ColAmd
  | Natural

  type info = {
    mutable ordering : ordering option;

    mutable reinit : int -> int option -> unit;
    mutable set_ordering : ordering -> unit;
  }

  let info = {
    ordering     = None;
    reinit       = (fun _ _ -> ());
    set_ordering = (fun _ -> ());
  }
end

module Superlumt = struct

  (* Must correspond with lsolver_ml.h:lsolver_superlumt_ordering_tag *)
  type ordering =
    Natural
  | MinDegreeProd
  | MinDegreeSum
  | ColAmd

  type info = {
    mutable ordering     : ordering option;
    mutable set_ordering : ordering -> unit;
  }

  let info = {
    ordering     = None;
    set_ordering = (fun _ -> ());
  }
end

module Direct = struct

  type solver =
    | Dense
    | LapackDense
    | Band
    | LapackBand
    | Klu of Klu.info
    | Superlumt of Superlumt.info

  type cptr

  type ('m, 'nd, 'nk) t = {
    rawptr : cptr;
    compat : solver;
  }
end

module Iterative = struct

  type info = {
    mutable maxl             : int option;
    mutable gs_type          : Spils.gramschmidt_type option;
    mutable max_restarts     : int option;

    mutable set_maxl         : int -> unit;
    mutable set_gs_type      : Spils.gramschmidt_type -> unit;
    mutable set_max_restarts : int -> unit;
  }

  let info = {
    maxl             = None;
    gs_type          = None;
    max_restarts     = None;

    set_maxl         = (fun _ -> ());
    set_gs_type      = (fun _ -> ());
    set_max_restarts = (fun _ -> ());
  }

  (* Must correspond with lsolver_ml.h:lsolver_iterative_solver_tag *)
  type solver =
    | Spbcgs
    | Spfgmr
    | Spgmr
    | Sptfqmr
    | Pcg

  type cptr

  type ('nd, 'nk, 'f) t = {
    rawptr : cptr;
    solver : solver;
    compat : info;
  }

  (* Must correspond with lsolver_ml.h:lsolver_gramschmidt_type_tag *)
  type gramschmidt_type = Spils.gramschmidt_type =
    | ModifiedGS
    | ClassicalGS

  (* Must correspond with lsolver_ml.h:lsolver_preconditioning_type_tag *)
  type preconditioning_type = Spils.preconditioning_type =
    | PrecNone
    | PrecLeft
    | PrecRight
    | PrecBoth

  external c_set_prec_type : cptr -> solver -> preconditioning_type -> unit
    = "ml_lsolvers_set_prec_type"

end

