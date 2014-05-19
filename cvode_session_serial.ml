(***********************************************************************)
(*                                                                     *)
(*               OCaml interface to (serial) Sundials                  *)
(*                                                                     *)
(*  Timothy Bourke (Inria), Jun Inoue (Inria), and Marc Pouzet (LIENS) *)
(*                                                                     *)
(*  Copyright 2014 Institut National de Recherche en Informatique et   *)
(*  en Automatique.  All rights reserved.  This file is distributed    *)
(*  under a BSD 2-Clause License, refer to the file LICENSE.           *)
(*                                                                     *)
(***********************************************************************)

include Cvode

type nvec = Sundials.Carray.t
type val_array = Sundials.Carray.t
type der_array = Sundials.Carray.t

type root_array = Sundials.Roots.t
type root_val_array = Sundials.Roots.val_array

type single_tmp = nvec
type triple_tmp = val_array * val_array * val_array

type 't jacobian_arg =
  {
    jac_t   : float;
    jac_y   : val_array;
    jac_fy  : val_array;
    jac_tmp : 't
  }

type prec_solve_arg =
  {
    rhs   : val_array;
    gamma : float;
    delta : float;
    left  : bool;
  }

(* Note this definition differs from the one in cvode_nvector, so the tag
   values are different.  *)
type linear_solver =
  | Dense of dense_jac_fn option
  | LapackDense of dense_jac_fn option
  | Band of bandrange * band_jac_fn option
  | LapackBand of bandrange * band_jac_fn option
  | Diag
  | Spgmr of spils_params * spils_callbacks
  | Spbcg of spils_params * spils_callbacks
  | Sptfqmr of spils_params * spils_callbacks
  | BandedSpgmr of spils_params * bandrange
  | BandedSpbcg of spils_params * bandrange
  | BandedSptfqmr of spils_params * bandrange

and dense_jac_fn = triple_tmp jacobian_arg -> Dls.DenseMatrix.t -> unit

and band_jac_fn = triple_tmp jacobian_arg -> int -> int -> Dls.BandMatrix.t -> unit

and bandrange = { mupper : int; mlower : int; }

and spils_params = { prec_type : Spils.preconditioning_type;
                     maxl : int option; }

and spils_callbacks =
  {
    prec_solve_fn : (single_tmp jacobian_arg -> prec_solve_arg -> nvec
                     -> unit) option;
    prec_setup_fn : (triple_tmp jacobian_arg -> bool -> float -> bool) option;
    jac_times_vec_fn : (single_tmp jacobian_arg -> val_array -> val_array
                        -> unit) option;
  }

type iter =
  | Newton of linear_solver
  | Functional

type cvode_mem
type cvode_file
type c_weak_ref

type session = {
        cvode      : cvode_mem;
        backref    : c_weak_ref;
        neqs       : int;
        nroots     : int;
        err_file   : cvode_file;

        mutable exn_temp   : exn option;

        mutable rhsfn      : float -> val_array -> der_array -> unit;
        mutable rootsfn    : float -> val_array -> root_val_array -> unit;
        mutable errh       : Sundials.error_details -> unit;
        mutable errw       : val_array -> nvec -> unit;
        mutable jacfn      : triple_tmp jacobian_arg -> Dls.DenseMatrix.t -> unit;
        mutable bandjacfn  : triple_tmp jacobian_arg -> int -> int
                               -> Dls.BandMatrix.t -> unit;
        mutable presetupfn : triple_tmp jacobian_arg -> bool -> float -> bool;
        mutable presolvefn : single_tmp jacobian_arg -> prec_solve_arg -> nvec
                               -> unit;
        mutable jactimesfn : single_tmp jacobian_arg -> val_array -> val_array
                               -> unit;

        mutable sensext    : Obj.t option (* Used by CVODES *)
      }

