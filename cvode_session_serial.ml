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

(* basic cvode types *)

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

type bandrange = { mupper : int; mlower : int; }

type spils_params = { maxl : int option;
                      prec_type : Spils.preconditioning_type; }

type dense_jac_fn = triple_tmp jacobian_arg -> Dls.DenseMatrix.t -> unit

type band_jac_fn = triple_tmp jacobian_arg -> int -> int -> Dls.BandMatrix.t -> unit

type spils_callbacks =
  {
    prec_solve_fn : (single_tmp jacobian_arg -> prec_solve_arg -> nvec
                     -> unit) option;
    prec_setup_fn : (triple_tmp jacobian_arg -> bool -> float -> bool) option;
    jac_times_vec_fn : (single_tmp jacobian_arg -> val_array -> val_array
                        -> unit) option;
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

type iter =
  | Newton of linear_solver
  | Functional

(* basic cvodes types *)

type quadrhsfn = float -> val_array -> der_array -> unit

type sensrhsfn =
    AllAtOnce of (float -> val_array -> der_array
                        -> val_array array -> der_array array
                        -> nvec -> nvec -> unit)
  | OneByOne of (float -> val_array -> der_array -> int
                       -> val_array -> der_array -> nvec -> nvec -> unit)

type quadsensrhsfn =
   float -> val_array -> val_array -> der_array -> val_array array
         -> nvec -> nvec -> unit

type brhsfn =
        BackBasic of (float -> val_array -> val_array -> der_array -> unit)
      | BackWithSens of (float -> val_array -> val_array array
                               -> val_array -> der_array -> unit)

type bquadrhsfn =
        QuadBasic of (float -> val_array -> val_array -> der_array -> unit)
      | QuadWithSens of (float -> val_array -> val_array array
                               -> val_array -> der_array -> unit)

type 't bjacobian_arg =
  {
    jac_t   : float;
    jac_y   : val_array;
    jac_yb  : val_array;
    jac_fyb : val_array;
    jac_tmp : 't
  }

type bprec_solve_arg =
  {
    rvecB  : val_array;
    gammaB : float;
    deltaB : float;
  }

type bdense_jac_fn =
      triple_tmp bjacobian_arg -> Dls.DenseMatrix.t -> unit

 and bband_jac_fn =
      bandrange -> triple_tmp bjacobian_arg -> Dls.DenseMatrix.t -> unit

(* the session type *)

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

        mutable sensext    : sensext (* Used by CVODES *)
      }

and sensext =
      NoSensExt
    | FwdSensExt of fsensext
    | BwdSensExt of bsensext

and fsensext = {
    (* Quadrature *)
    mutable quadrhsfn       : quadrhsfn;

    (* Forward *)
    mutable num_sensitivies : int;
    mutable sensarray1      : val_array array;
    mutable sensarray2      : der_array array;
    mutable senspvals       : Sundials.real_array option;
                            (* keep a reference to prevent garbage collection *)

    mutable sensrhsfn       : (float -> val_array -> der_array -> val_array array
                               -> der_array array -> nvec -> nvec -> unit);
    mutable sensrhsfn1      : (float -> val_array -> der_array -> int -> val_array 
                               -> der_array -> nvec -> nvec -> unit);
    mutable quadsensrhsfn   : quadsensrhsfn;
  }

and bsensext = {
    (* Adjoint *)
    parent                : session ;
    which                 : int;

    bsensarray            : val_array;

    mutable brhsfn        : (float -> val_array -> val_array
                                   -> der_array -> unit);
    mutable brhsfn1       : (float -> val_array -> val_array array
                                   -> val_array -> der_array -> unit);
    mutable bquadrhsfn    : (float -> val_array -> val_array
                                   -> der_array -> unit);
    mutable bquadrhsfn1   : (float -> val_array -> val_array array
                                   -> val_array -> der_array -> unit);

    mutable bpresetupfn : triple_tmp bjacobian_arg -> bool -> float -> bool;
    mutable bpresolvefn : single_tmp bjacobian_arg
                            -> bprec_solve_arg -> nvec -> unit;
    mutable bjactimesfn : single_tmp bjacobian_arg -> nvec -> nvec -> unit;

    mutable bjacfn      : bdense_jac_fn;
    mutable bbandjacfn  : bband_jac_fn;
  }

let shouldn't_be_called fcn =
  failwith ("internal error in sundials: " ^ fcn ^ " is called")
let dummy_dense_jac _ _ = shouldn't_be_called "dummy_dense_jac"
let dummy_band_jac _ _ _ _ = shouldn't_be_called "dummy_band_jac"
let dummy_prec_setup _ _ _ = shouldn't_be_called "dummy_prec_setup"
let dummy_prec_solve _ _ _ = shouldn't_be_called "dummy_prec_solve"
let dummy_jac_times_vec _ _ _ = shouldn't_be_called "dummy_jac_times_vec"
let dummy_bprec_setup _ _ _ = shouldn't_be_called "dummy_prec_setup"
let dummy_bprec_solve _ _ _ = shouldn't_be_called "dummy_prec_solve"
let dummy_bjac_times_vec _ _ _ = shouldn't_be_called "dummy_jac_times_vec"
let dummy_bdense_jac _ _ = shouldn't_be_called "dummy_dense_jac"
let dummy_bband_jac _ _ _ _ = shouldn't_be_called "dummy_band_jac"

